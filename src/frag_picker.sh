#!/bin/bash

#BASE=/home/h3nnn4n/rosetta/rosetta_bin_linux_2017.26.59567_bundle
#BASE=/mnt/c/Users/Nilcimar/Desktop/protein-prediction-framework/rosetta_bin_linux_2019.35.60890_bundle
BASE=/home/nilcimar/workspace/protein-prediction-framework/rosetta_bin_linux_2019.35.60890_bundle
VALL=$BASE/tools/fragment_tools/vall.jul19.2011.gz
DATABASE=$BASE/main/database
BIN=$BASE/main/source/bin/fragment_picker.static.linuxgccrelease

NAME=$1
FASTA=$NAME.fasta
NATIVE=$NAME.pdb
PSIPRED=$NAME.psipred.ss2
SPIDER=$NAME.spider.ss2
RAPTORX=$NAME.raptorx.ss2
PSSM=$NAME.pssm
#MUFOLD=$NAME.mufold.ss2

OUT=output

if [ ! -d $OUT ]
then
    mkdir $OUT
fi

# utilizar este comando para gerar os fragmentos com o protocolo best
# substituir $RAPTORX pela variavel que representa o preditor de estrutura secundaria desejado
#$BIN -database $DATABASE -in::file::vall $VALL -in::file::fasta $FASTA -in::file::s $NATIVE -in:file:pssm $PSSM -frags::ss_pred $RAPTORX predA -frags::denied_pdb $NAME.homolog_vall -frags::scoring::config "../../src/scoring.wghts" -frags::bounded_protocol -frags::frag_sizes 9 3 -frags::n_candidates 200 -frags::n_frags 200 -out::file::frag_prefix $OUT/$NAME -frags::describe_fragments $OUT/$NAME.fsc

# utilizar este comando para gerar os fragmentos com o protocolo quota
$BIN -database $DATABASE -in::file::vall $VALL -in::file::fasta $FASTA -in::file::s $NATIVE -in:file:pssm $PSSM -frags::ss_pred $RAPTORX predA $SPIDER predB $PSIPRED predC -frags::denied_pdb $NAME.homolog_vall -out::file::frag_prefix $OUT/$NAME -frags::describe_fragments $OUT/$NAME.fsc -frags::picking::quota_config_file "../../src/quota.def" -frags::scoring::config "../../src/quota_protocol.wghts" -frags::frag_sizes 9 3 -frags::n_candidates 1000 -frags::n_frags 200 

#-ignore_unrecognized_res

#-frags::denied_pdb $NAME.homolog_vall