#!/bin/bash

PLINKOUT=/mnt/data1/lawless/spss/plink
COUNT=$PLINKOUT/count_varID_length/

# write just the variant ID
for file in $COUNT/*varID_over79; do
    egrep -o '^[^=]*' $file > $file.original_longID
done 

# replace anything after second : delim
for file in $COUNT/*varID_over79; do
    sed 's/:[^:]*/:ID/2g' $file > $file.ID
done


regexp=$(printf %s "$old" | sed 's:[$*./\[^]:\\&:g')
replacement=$(printf %s "$new" | sed 's:[\&/]:\\&:g')
sed -e "s/$regexp/$replacement/g"

for each file in the filelist
    for each line in list
        sed replcace

old = $COUNT/$file.original_longID
new = $COUNT/$file.ID


cust_func(){
	while read line; do
	  awk '{ print $2 " = " length($2) }' $PLINKOUT/$line.bim > $PLINKOUT/count_varID_length/$line.count
	done < $FILELIST/chr$i\_files #| while read line; do
}

#for i in 23 24; do
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
cust_func &
done
wait

