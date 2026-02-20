#!/bin/bash

DATADIR="../data/Wheat_Relative_History_Data_Glemin_et_al/OneCopyGenes-trimmed"

mkdir ../results/
mkdir ../results/RAxML/

for file in "$DATADIR"/*; do
    raxml-ng --msa $file --model GTR+G4
    mv ${file}**raxml** ../results/RAxML/
done