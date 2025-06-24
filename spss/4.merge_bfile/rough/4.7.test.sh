#!/bin/bash

PLINKOUT=/mnt/data1/lawless/spss/plink
COUNT=$PLINKOUT/count_varID_length/

# while read line; do
# Grep -w $PLINKOUT/SPSS.QC_chr9.pos95688099-100684757.impute2.bim -e $line
# Done < $COUNT/Var.originalID


filename=Var.original
bimfile=$PLINKOUT/$filename.bim
replacement=$PLINKOUT/$filename.replacedID.bim
match=$filename.ID

# replace anything after second : delim
while read line; do
# if
#     grep -w $bimfile -e $line; then
     sed 's/$line/newID/g' $bimfile > $replacement
# fi
done < $COUNT/$match

