import subprocess
import requests
import shutil
import yaml
import gzip
import Bio.PDB
import sys
import os


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

    if os.path.exists(target[:4].lower()):
        print("Target already %s exists on %s, skipping..." % (target[:4].lower(), config['output_path']))
        os.chdir(original_path)
        return

    print("Summoning %s" % target)
    pdbid, chain = get_pdb_code_and_chain(target)
    target = pdbid.lower()

    os.makedirs(target)
    os.chdir(target)

    fetch_native_pdb(target, pdbid)
    fetch_fasta(target, pdbid, chain)

    if config['do_psspred']:
        os.makedirs('psspred')
        run_psspred(config, target)

    if config['do_psipred']:
        run_psipred(config, target)

    if config['do_fragpicking']:
        os.makedirs('output')
        frag_path = config['frag_path']
        subprocess.call([frag_path, target])

    os.chdir(original_path)


def get_pdb_code_and_chain(target):
    pdbid = target[0:4].upper()
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
    r = requests.get(url)

    with open(target + '.pdb.gz', 'wb') as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)

    with gzip.open(target + '.pdb.gz', 'rb') as f_in, open(target + '.pdb', 'wb') as f_out:
        f_out.write(f_in.read())


def fetch_fasta(target, pdbid, chain):
    url = "http://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=%s&compressionType=uncompressed" % pdbid
    r = requests.get(url)

    with open(target + '.fasta', 'wb') as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)

    chain_count = count_nchain_in_fasta(target)

    if chain_count > 1:
        if chain is None:
            raise Exception('Target %s has %d chains but not was specified' % (pdbid, chain_count))
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

    os.chdir('psspred')

    subprocess.call([pss_path, '../' + target + '.fasta'])

    os.chdir("../")

    with open('psspred/seq.dat.ss', 'rt') as f_in, open(target + '.psipred.ss2', 'wt') as f_out:
        c = 0
        for line in f_in.readlines():
            if c == 0:
                c += 1
                f_out.write('# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n')
            else:
                f_out.write(line)


def run_psipred(config, target):
    psipred_path = config['psipred_path']
    original_path = os.getcwd()
    os.chdir(psipred_path)

    subprocess.call([
        os.path.join(psipred_path, config['psipred_bin']),
        os.path.join(original_path, (target + '.fasta'))
    ])

    os.chdir(original_path)

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
