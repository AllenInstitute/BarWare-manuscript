#!/bin/bash

# profile processing time, peak memory usage of Cellrange Count
printf "script %s is being run by %s with PID %s\n" $0 $USER $$

# read list of wells
MANIFESTS="/mnt/disks/barcode-tender-manuscript/barcode_manifests"
LIBRARIES="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/cellranger/library_manifests"
WELLS=$( ls $MANIFESTS | grep .tsv.gz )
TAGLIST="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/cellranger/cell_hashing_feature-ref.csv"
REFERENCE="/shared/apps/refdata-gex-GRCh38-2020-A"

# initialize variables to NULL
NAME=""
LIB_SHEET=""

# call command for each library
for W in $WELLS; do
    NAME=$( echo $W | sed s'/_barcodes.*//' )
    LIB_SHEET=$( ls -d $LIBRARIES/* | grep $NAME )

    # set command
    COMMAND="/shared/apps/cellranger-4.0.0/cellranger count --nosecondary --nopreflight --disable-ui --transcriptome=${REFERENCE} --id=${NAME} --libraries=${LIB_SHEET} --feature-ref=${TAGLIST} --expect-cells=40000"

    echo "Runing Cellranger Count on ${NAME}"
    printf "\n Command: \n\n$COMMAND\n"

    # run and time command
    $( command time -v -o "/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/cellranger/CellrangerCount_${NAME}_profile.txt" $COMMAND )

    echo "counting complete for $NAME"
    echo

done

printf "Cellranger profiling complete\n"
