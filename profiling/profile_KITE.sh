#!/bin/bash

# profile processing time, peak memory usage of KITE
printf "script %s is being run by %s with PID %s\n" $0 $USER $$

# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

# GitHub link:
# https://github.com/pachterlab/kite

# -- previous to running this pipeline the taglist was converted from CSV format to FASTA and t2g files as instructed in the GitHub README using the command below:

# "python3 /shared/apps/kite/featuremap/featuremap.py KITE_cell_hashing_taglist.csv --t2g $PWD/Cell_Hash_Mismatch.t2g --fa $PWD/Cell_Hash_Mismatch.fa --header"

# TAGLIST: "/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/kite_cell_hashing_taglist.csv"
# FASTA file: "/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/Cell_Hash_Mismatch.fa"
# t2g file: "/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/Cell_Hash_Mismatch.t2g"

# -- next, a kallisto index was build using the following command:
# "/shared/apps/kallisto/kallisto index -i Cell_Hash_Mismatch.idx -k 15 ./Cell_Hash_Mismatch.fa"

# Kallisto Index: "/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/Cell_Hash_Mismatch.idx"

# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------


# read list of wells
MANIFESTS="/mnt/disks/barcode-tender-manuscript/barcode_manifests"
WELLS=$( ls $MANIFESTS | grep .tsv.gz )
WHITELIST="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/3M-february-2018.txt"

# initialize variables to NULL
NAME=""

# call command for each library
for W in $WELLS; do
    NAME=$( echo $W | sed s'/_barcodes.*//' )

    # set command
    COMMAND="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/profile_KITE_helper.sh ${NAME}"

    echo "Running KITE pipeline for ${NAME}"

    # run and time command
    $( command time -v -o "/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/kite_${NAME}_profile.txt" $COMMAND )

    echo "KITE profiling complete for $NAME"
    echo

done

printf "KITE profiling complete\n"
