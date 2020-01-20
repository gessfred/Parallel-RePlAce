#!/bin/bash


echo "Cleaning output dir..."
rm -r ../output/etc
mkdir ../output/etc

for benchmark in $@
do
        for T in 1 2 4 8 12
        do
                ./bookshelf.sh $T $benchmark "baseline.csv"
        done
done
