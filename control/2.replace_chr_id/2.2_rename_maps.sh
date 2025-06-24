#!/bin/bash

cd /work/gr-fe/lawless/spss/gwas/shcs

for file in *.chr.map; do
  mv "$file" "${file%.chr.map}.map"
done

