## download all exons
```
GENOME=hg19
curl  -s "http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/database/knownGene.txt.gz" | gunzip -c |\
 awk '{n=int($8); split($9,S,/,/);split($10,E,/,/); for(i=1;i<=n;++i) {printf("%s\t%s\t%s\t%s\t%s\n",$2,S[i],E[i],$3,$1);} }' >\
 exons.${GENOME}.bed
```
