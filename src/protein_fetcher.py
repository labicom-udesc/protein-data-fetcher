import subprocess
import requests
import shutil
import yaml
import gzip
#import Bio.PDB
import sys
import os
import re


class Selector:
    def __init__(self, chain):
        self.chain = chain.upper()

    def accept_chain(self, chain):
        return (chain.get_id().upper() == self.chain)

    def accept_model(self, _):
        return True

    def accept_residue(self, _):
        return True

    def accept_atom(self, _):
        return True


def fetch(config, target):
    if not os.path.exists(config['output_path']):
        os.makedirs(config['output_path'])

    original_path = os.getcwd()
    os.chdir(config['output_path'])

    '''if os.path.exists(target[:4].lower()):
        print("Target already %s exists on %s, skipping..." % (target[:4].lower(), config['output_path']))
        os.chdir(original_path)
        return'''

    print("Summoning %s" % target)
    pdbid, chain = get_pdb_code_and_chain(target)
    target = pdbid.lower()

    #os.makedirs(target)
    os.chdir(target)

    #fetch_native_pdb(target, pdbid)
    #fetch_fasta(target, pdbid, chain)

    if config['do_sspro8']:
        os.makedirs('sspro8')
        #run_psspred(config, target)

    if config['do_raptorx']:
        #os.makedirs('raptorx')
        run_raptorx(config, target)

    if config['do_spider']:
        #os.makedirs('spider')
        run_spider(config, target)

    if config['do_mufold']:
        #os.makedirs('spider')
        run_mufold(config, target)

    if config['do_psspred']:
        os.makedirs('psspred')
        run_psspred(config, target)

    if config['do_psipred']:
        os.makedirs('psipred')
        run_psipred(config, target)

    if config['do_fragpicking']:
        os.makedirs('output')
        frag_path = config['frag_path']
        subprocess.call([frag_path, target])

    os.chdir(original_path)


def get_pdb_code_and_chain(target):
    #pdbid = target[0:4].upper()
    pdbid = target[0:5].upper()
    chain = None

    if '_' in target:
        chain = target.split('_')[1]

    print("using pdbid: %4s" % pdbid, end='')
    if chain is not None:
        print(' and chain: %s' % chain, end='')
    print()

    return pdbid, chain


def fetch_native_pdb(target, pdbid):
    url = "http://files.rcsb.org/download/%s.pdb.gz" % pdbid
    #url = "http://files.rcsb.org/download/4hhb.pdb.gz"
    r = requests.get(url)

    with open(target + '.pdb.gz', 'wb') as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)

    with gzip.open(target + '.pdb.gz', 'rb') as f_in, open(target + '.pdb', 'wb') as f_out:
        f_out.write(f_in.read())


def fetch_fasta(target, pdbid, chain):
    url = "https://www.rcsb.org/fasta/entry/%s" % pdbid
    #url = "http://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=%s&compressionType=uncompressed" % pdbid
    #url = "http://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=4hhb&compressionType=uncompressed"
    r = requests.get(url)

    with open(target + '.fasta', 'wb') as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)

    chain_count = count_nchain_in_fasta(target)

    if chain_count > 1:
        if chain is None:
            raise Exception('Target %s has %d chains but not was specified' %(pdbid, chain_count))
        extract_chain_from_fasta(target, pdbid, chain)
        extract_chain_from_pdb(target, pdbid, chain)


def extract_chain_from_fasta(target, pdbid, target_chain):
    lines = ''
    with open(target + '.fasta', 'rt') as fd:
        for line in fd.readlines():
            lines = lines + line

    for line in lines.split('>'):
        if len(line) > 5:
            chain = line[5]
            if chain == target_chain:
                with open(target + '.fasta', 'wt') as fd:
                    fd.write('>' + line)
                return

    raise Exception('Chain %s was not found on %s' % (chain, pdbid))


def extract_chain_from_pdb(target, pdbid, chain):
    parser = Bio.PDB.PDBParser()
    pdbio = Bio.PDB.PDBIO()
    chain_pdb = parser.get_structure(pdbid, pdbid.lower() + '.pdb')
    pdbio.set_structure(chain_pdb)
    pdbio.save(target + '.pdb', select=Selector(chain))


def count_nchain_in_fasta(target):
    chain_count = 0
    with open(target + '.fasta', 'rt') as fd:
        for line in fd.readlines():
            if '>' in line:
                chain_count += 1

    return chain_count


def run_psspred(config, target):
    pss_path = config['pss_path']
    original_path = os.getcwd()

    os.chdir('psspred')

    subprocess.call([
        pss_path,
        os.path.join(original_path, (target + '.fasta'))
    ])

    #subprocess.call([pss_path, '../' + target + '.fasta'])

    os.chdir("../")

    with open('psspred/seq.dat.ss', 'rt') as f_in, open(target + '.psspred.ss2', 'wt') as f_out:
        c = 0
        for line in f_in.readlines():
            if c == 0:
                c += 1
                f_out.write('# PSIPRED VFORMAT (PSIPRED V4.0)\n')
            else:
                f_out.write(line)


def run_raptorx(config, target):
    #spider_path = config['spider_path']
    original_path = os.getcwd()

    os.chdir('raptorx')

    os.chdir("../")

    with open('raptorx/' + target + '.ss3', 'rt') as f_in, open(target + '.raptorx.ss2', 'wt') as f_out:
        c = 0
        t = 0
        for line in f_in.readlines():
            s = line
            s = s.replace("   ", ',')
            s = s.replace("  ", ',')
            s = s.replace(' ', ',')
            #print(s)
            entries = re.split(',', s.strip())
            #print(entries)
            if c == 0:
                f_out.write('# PSIPRED VFORMAT (PSIPRED V4.0)\n\n')
            elif c > 2:
                t = t+1
                ps = entries[2]
                ss = entries[3]
                ss_c = float(entries[6])
                ss_e = float(entries[5])
                ss_h = float(entries[4])
                #linha = t + ' ' + ps + ' ' + ss + '  ' + ss_c + '  ' + ss_h + '  ' + ss_e + '\n'
                f_out.write("%4d %c %c   %5.3f  %5.3f  %5.3f\n" % (t,ps,ss,ss_c,ss_h,ss_e))
            c += 1

def run_spider(config, target):
    #spider_path = config['spider_path']
    original_path = os.getcwd()

    os.chdir('spider')

    os.chdir("../")

    with open('spider/' + target + '.spd3', 'rt') as f_in, open(target + '.spider.ss2', 'wt') as f_out:
        c = 0
        t = 0
        for line in f_in.readlines():
            entries = re.split('\t', line.strip())
            #print(entries)
            if c == 0:
                c += 1
                f_out.write('# PSIPRED VFORMAT (PSIPRED V4.0)\n\n')
            else:
                t = t+1
                ps = entries[1]
                ss = entries[2]
                ss_c = float(entries[8])
                ss_e = float(entries[9])
                ss_h = float(entries[10])
                #linha = t + ' ' + ps + ' ' + ss + '  ' + ss_c + '  ' + ss_h + '  ' + ss_e + '\n'
                f_out.write("%4d %c %c   %5.3f  %5.3f  %5.3f\n" % (t,ps,ss,ss_c,ss_h,ss_e))


def run_mufold(config, target):
    original_path = os.getcwd()
    print(original_path)

    os.chdir('mufold')

    os.chdir("../")

    arq = open('mufold/query.fasta','r')
    fasta = arq.readlines()[1].strip()

    arq2 = open('mufold/query.ss','r')
    seq = arq2.readlines()
    ss3 = seq[1].strip()
    ss8 = seq[3].strip()
    print(target)

    with open('mufold/query.fasta.q3prob', 'rt') as f_in, open(target + '.mufold.ss2', 'wt') as f_out:
        c = 0
        t = 0
        for line in f_in.readlines():
            entries = re.split(',', line.strip())
            print(entries)
            if c == 0:
                c += 1
                f_out.write('# PSIPRED VFORMAT (PSIPRED V4.0)\n\n')
            else:
                t = t+1
                ps = fasta[t-1]
                ss = ss3[t-1]
                ss_c = float(entries[2])
                ss_e = float(entries[1])
                ss_h = float(entries[0])
                f_out.write("%4d %c %c   %5.3f  %5.3f  %5.3f\n" % (t,ps,ss,ss_c,ss_h,ss_e))


def run_psipred(config, target):
    psipred_path = config['psipred_path']
    original_path = os.getcwd()
    os.chdir(psipred_path)

    subprocess.call([
        os.path.join(psipred_path, config['psipred_bin']),
        os.path.join(original_path, (target + '.fasta'))
    ])

    os.chdir(original_path)
    print('\n\n\nEXECUTANDO\n\n\n')

    src_file = os.path.join(psipred_path, target + '.ss2')
    dest_file = os.path.join(original_path, target + '.psipred.ss2')
    shutil.copyfile(src_file, dest_file)


def config_loader():
    config_path = './protein_fetcher_config.yaml'

    with open(config_path, 'rt') as f:
        config = yaml.load(f)

    return config


if __name__ == "__main__":
    config = config_loader()

    args = sys.argv[1:]
    if len(args) < 1:
        print("Missing arguments")

    for target in args:
        fetch(config, target)
