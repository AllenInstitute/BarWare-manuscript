import pandas as pd
import os
from glob import glob
import csv

os.chdir("/mnt/disks/barcode-tender-manuscript/cell-hashing-pipeline/analysis")

# paths to hto_parser results for each calling method
barcounter_calls="/mnt/disks/barcode-tender-manuscript/cell-hashing-pipeline/hto_parser/barcounter"
CITE_calls="/mnt/disks/barcode-tender-manuscript/cell-hashing-pipeline/hto_parser/CITE-seq_Count"

# output names
names = ["Pool-16","Pool-24","Pool-32","Pool-48","Pool-64","Pool-80"]

# store categories in list of dataframes
barcounter_cats = []
for n in names:
    barcounter_cats.append(pd.read_csv(f'{barcounter_calls}/{n}_parsing_results/{n}_hto_category_table.csv.gz'))

cite_cats = []
for n in names:
    cite_cats.append(pd.read_csv(f'{CITE_calls}/{n}_parsing_results/{n}_hto_category_table.csv.gz'))

# format header for output CSV file
bar_cite_overlap = []
total_counts = []

# find overlap with singlet, doublet calls --------------------

# total barcounter singlets acrross all wells
barcounter_s = 0
for i in range(len(names)):
    barcounter_s += len(barcounter_cats[i][barcounter_cats[i]['hto_category'] == 'singlet']['cell_barcode'])

# total singlets shared between both methods
total_s_overlap = 0
for i in range(len(names)):
    s1 = tuple(barcounter_cats[i][barcounter_cats[i]['hto_category'] == 'singlet']['cell_barcode'])
    s2 = set(cite_cats[i][cite_cats[i]['hto_category'] == 'singlet']['cell_barcode'])
    for bc in s1:
        if bc in s2:
            total_s_overlap += 1

# CITE-seq_Count_Correction
temp = ['singlet']
for i in range(len(names)):
    overlap = 0
    s1 = tuple(barcounter_cats[i][barcounter_cats[i]['hto_category'] == 'singlet']['cell_barcode'])
    s2 = set(cite_cats[i][cite_cats[i]['hto_category'] == 'singlet']['cell_barcode'])
    total_counts.append((len(s1),len(s2)))
    for bc in s1:
        if bc in s2:
            overlap += 1
    temp.append(overlap / len(s1))
bar_cite_overlap.append(temp)

temp = ['doublet']
for i in range(len(names)):
    overlap = 0
    s1 = tuple(barcounter_cats[i][barcounter_cats[i]['hto_category'] == 'doublet']['cell_barcode'])
    s2 = set(cite_cats[i][cite_cats[i]['hto_category'] == 'doublet']['cell_barcode'])
    total_counts.append((len(s1),len(s2)))
    for bc in s1:
        if bc in s2:
            overlap += 1
    temp.append(overlap / len(s1))
bar_cite_overlap.append(temp)

temp = ['multiplet']
for i in range(len(names)):
    overlap = 0
    s1 = tuple(barcounter_cats[i][barcounter_cats[i]['hto_category'] == 'multiplet']['cell_barcode'])
    s2 = set(cite_cats[i][cite_cats[i]['hto_category'] == 'multiplet']['cell_barcode'])
    total_counts.append((len(s1),len(s2)))
    for bc in s1:
        if bc in s2:
            overlap += 1
    temp.append(overlap / len(s1))
bar_cite_overlap.append(temp)

temp = ['no_hash']
for i in range(len(names)):
    overlap = 0
    s1 = tuple(barcounter_cats[i][barcounter_cats[i]['hto_category'] == 'no_hash']['cell_barcode'])
    s2 = set(cite_cats[i][cite_cats[i]['hto_category'] == 'no_hash']['cell_barcode'])
    total_counts.append((len(s1),len(s2)))
    for bc in s1:
        if bc in s2:
            overlap += 1
    temp.append(overlap / len(s1))
bar_cite_overlap.append(temp)


# compare population identity between BarCounter and CITE-seq Count w/ UMI correction for each mixed well --------------------

temp = ['singlet_pops']
for i in range(len(names)):
    assoc1 = dict()
    assoc2 = dict()
    overlap = 0
    total = 0
    bc1 = list(barcounter_cats[i][barcounter_cats[i]['hto_category'] == 'singlet']['cell_barcode'])
    pop1 = list(barcounter_cats[i][barcounter_cats[i]['hto_category'] == 'singlet']['pbmc_sample_id'])
    bc2 = list(cite_cats[i][cite_cats[i]['hto_category'] == 'singlet']['cell_barcode'])
    pop2 = list(cite_cats[i][cite_cats[i]['hto_category'] == 'singlet']['pbmc_sample_id'])

    # translate population calls into dict
    for i in range(len(bc1)):
        assoc1[bc1[i]] = pop1[i]
    for i in range(len(bc2)):
        assoc2[bc2[i]] = pop2[i]

    # cross-check each barcode between population calls
    for k,v in assoc1.items():
        if k in assoc2.keys():
            total += 1
            if assoc2[k] == v:
                overlap += 1

    temp.append(overlap / total)

bar_cite_overlap.append(temp)

# output results to CSV file
with open('barcode_call_overlaps.csv','w') as fw:
    csv_writer = csv.writer(fw)
    csv_writer.writerow(['comparison'] + names)
    for row in bar_cite_overlap:
        csv_writer.writerow(row)


# look into disagreements between BarCounter and CITE-seq Count w/ UMI correction --------------------

# store count matrices in list of dataframes
barcounter_counts = []
for n in names:
    barcounter_counts.append(pd.read_csv(f'{barcounter_calls}/{n}_parsing_results/{n}_hto_count_matrix.csv.gz'))

cite_counts = []
for n in names:
    cite_counts.append(pd.read_csv(f'{CITE_calls}/{n}_parsing_results/{n}_hto_count_matrix.csv.gz'))

# get list of all BarCounter singlets called as doublets by CITE-seq Count (with UMI correction)
swap_header = ['Well','Barcode','BarCounter_1st','BarCounter_2nd','CITE_1st','CITE_2nd','Change_1st','Change_2nd','Prop_Change_1st','Prop_change_2nd','Ratio_1st:2nd']
swaps = []
for i in range(len(names)):
    s1 = tuple(barcounter_cats[i][barcounter_cats[i]['hto_category'] == 'singlet']['cell_barcode'])
    s2 = set(cite_cats[i][cite_cats[i]['hto_category'] == 'doublet']['cell_barcode'])
    for bc in s1:
        if bc in s2:
            # find difference in ration between highest and second highest counts between methods
            bc_list = list(barcounter_counts[i][bc])
            cite_list = list(cite_counts[i][bc])

            bc_first = bc_list.pop(bc_list.index(max(bc_list)))
            bc_sec = bc_list.pop(bc_list.index(max(bc_list)))
            cite_first = cite_list.pop(cite_list.index(max(cite_list)))
            cite_sec = cite_list.pop(cite_list.index(max(cite_list)))

            swaps.append([names[i], bc, bc_first, bc_sec, cite_first, cite_sec, (bc_first - cite_first), (bc_sec - cite_sec),(bc_first - cite_first)/bc_first, (bc_sec - cite_sec)/bc_sec, (bc_first / bc_sec)])

# output results to CSV file
with open('doublet_swap_counts.csv','w') as fw:
    csv_writer = csv.writer(fw)
    csv_writer.writerow(swap_header)
    for row in swaps:
        csv_writer.writerow(row)