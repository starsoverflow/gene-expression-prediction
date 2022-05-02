from os import listdir
from os.path import join, isfile

import torch
import collections
import pdb
import torch.utils.data
import csv
import json
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms, utils
import math
import pandas as pd
import os
import numpy as np

DATA_DIR = "./data/"
OUTPUT_DIR = "./processed_data/"

def clean_rna_data():
    esc_rna = pd.read_csv(DATA_DIR + "GSE109262_ESC_RNA_counts.txt", sep='\t')
    cols = esc_rna.columns

    valid_genes = []
    for idx, row in esc_rna.iterrows():
        gene = row["ESC_ens_id"]
        cnt = 0
        for col in cols:
            if col != "ESC_ens_id":
                cnt += row[col]
        if cnt > 1000:
            valid_genes.append(gene)

    esc_rna_out = esc_rna.loc[esc_rna["ESC_ens_id"].isin(valid_genes)]
    print(len(valid_genes))
    print(esc_rna_out.shape[0])
    esc_rna_out.to_csv(DATA_DIR + "RNA_counts_cleaned.txt", sep='\t', index=False)


def load_rna_counts():
    esc_rna = pd.read_csv(DATA_DIR + "RNA_counts_cleaned.txt", sep='\t')
    esc_rna_counts = {}

    cols = esc_rna.columns

    for idx, row in esc_rna.iterrows():
        gene = row["ESC_ens_id"]
        if gene not in esc_rna_counts:
            esc_rna_counts[gene] = {}
        for col in cols:
            if col != "ESC_ens_id":
                esc_rna_counts[gene][col] = row[col]

    return esc_rna_counts


def get_expr_data():
    cnt_dict = {}
    esc_rna_counts = load_rna_counts()
    max_in_gene = {}
    expressed_count = {}
    for gene in esc_rna_counts:
        for col in esc_rna_counts[gene]:
            if col not in cnt_dict:
                cnt_dict[col] = {}
            cnt_dict[col][gene] = esc_rna_counts[gene][col]
            if gene not in max_in_gene:
                max_in_gene[gene] = 0
            max_in_gene[gene] = max(max_in_gene[gene], esc_rna_counts[gene][col])

    for col in cnt_dict:
        file = open(OUTPUT_DIR + "Cell_"+col+".expr.csv", "w+")
        data = sorted(cnt_dict[col].items())
        for gene, cnt in data:
            if cnt >  max_in_gene[gene] / 4:
                value = 1
                if col not in expressed_count:
                    expressed_count[col] = 0
                expressed_count[col] += 1
            else:
                value = 0
            file.write(gene + "," + str(value) + "\n")

    cnts = np.array([v for k, v in expressed_count.items()])

    print(cnts.mean())
    print(cnts.max())
    print(cnts.min())

    # see this!
    # get top k cells
    arr = (sorted(expressed_count, key=expressed_count.get, reverse=True)[:10])
    print(arr)


def analyze_expr_data():
    cnt_dict = {}
    esc_rna_counts = load_rna_counts()
    for gene in esc_rna_counts:
        for col in esc_rna_counts[gene]:
            cnt_dict[col] = 0
        break

    for gene in esc_rna_counts:
        for col in esc_rna_counts[gene]:
            cnt_dict[col] += esc_rna_counts[gene][col]

    cnts = np.array([v for k, v in cnt_dict.items()])

    print(cnt_dict)
    print(cnts.mean())
    print(cnts.max())
    print(cnts.min())

    # get top k cells
    arr = (sorted(cnt_dict, key=cnt_dict.get, reverse=True)[:10])
    print(arr)

    files = [f for f in listdir(DATA_DIR+"/GSE109262_TSV") if f.endswith(".tsv")]

    for k in arr:
        for file in files:
            if file.find(k) != -1:
                print(k, cnt_dict[k], file)
                break


def clean_tss_data():
    esc_rna_counts = load_rna_counts()

    tss_data = pd.read_csv(DATA_DIR + "TSS.tsv", sep='\t')
    tss_data = tss_data.sort_values(['Gene stable ID', 'Transcription start site (TSS)', 'Chromosome/scaffold name'])
    previous_gene = ""
    data = {"gene": [], "tss": [], "chrom": []}
    for idx, row in tss_data.iterrows():
        gene = row["Gene stable ID"]
        if gene == previous_gene:
            continue
        if gene not in esc_rna_counts:
            continue
        tss = int(row["Transcription start site (TSS)"])
        chrom = row["Chromosome/scaffold name"]
        previous_gene = gene
        if chrom not in [str(x) for x in range(0, 20)]:
            continue
        data["tss"].append(tss)
        data["gene"].append(gene)
        data["chrom"].append(chrom)

    df = pd.DataFrame.from_dict(data)
    df = df.sort_values(["gene", "tss", "chrom"])
    df.to_csv(DATA_DIR + "TSS_clean.csv", sep='\t')


# clean_rna_data()
# get_expr_data()
# clean_tss_data()
analyze_expr_data()
# clean_rna_data()
