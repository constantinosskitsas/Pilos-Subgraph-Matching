#!/bin/bash
for i in $(seq 1 400)
do
  for j in 32 64 96 #128 256
  do
  	for t in dblp
	do
  		(timeout 600s ./SubgraphMatching.out -dataset $t -qsize $j -qnumber $i -qprop G -filter PL -alpha 125 -num 100000 -SF NOrd)
	done  
done
done