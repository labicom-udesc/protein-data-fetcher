#!/bin/bash

BASE=/home/h3nnn4n/rosetta/rosetta_bin_linux_2017.26.59567_bundle
VALL=$BASE/tools/fragment_tools/vall.jul19.2011.gz
DATABASE=$BASE/main/database
BIN=$BASE/main/source/bin/fragment_picker.static.linuxgccrelease

NAME=$1
FASTA=$NAME.fasta
NATIVE=$NAME.pdb
SSPRED=$NAME.psipred.ss2

OUT=output

if [ ! -d $OUT ]
then
    mkdir $OUT
fi

$BIN -database $DATABASE -in::file::vall $VALL -in::file::fasta $FASTA -in::file::s $NATIVE -frags::ss_pred $SSPRED predA -out::file::frag_prefix $OUT/$NAME -frags::describe_fragments $OUT/$NAME.fsc -frags::scoring::config "../../src/scoring.wghts" -ignore_unrecognized_res
