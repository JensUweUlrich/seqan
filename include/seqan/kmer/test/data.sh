#!/usr/bin/env bash

length=4194304
seed=20180406

for bin in 64 256 1024 8192
do
    binLength=$((length / bin))
    echo "Simulating $bin bins with reference length of $length and binLength of $binLength"
    mkdir -p $bin
    mkdir -p $bin/bins
    mkdir -p $bin/info
    mkdir -p $bin/reads
    cd $bin/bins
    # Simulate reference
    ../../mason_limited/mason_genome -l $length -o ref.fasta -s $seed &>/dev/null
    # Evenly distribute it over 64 bins
    pyfasta split -n $bin -k $binLength ref.fasta &>/dev/null
    # We do not need the reference anymore
    rm ref.fasta
    # Rename the bins
    for i in *.fasta; do mv $i bin_$(echo "$i" | cut -d "." -f2).fa; done
    # Simulate 16 haplotypes for each bin
    for i in *.fa
    do
        ../../mason_limited/mason_variator \
        -ir $i \
        -n 16 \
        -of $(basename $i .fa).fasta \
        -ov ../info/$(basename $i .fa).vcf &>/dev/null && rm $i
    done
    cd ../..
done
