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
from os import listdir
from os.path import isfile, join

DATA_DIR = "./processed_data/"

# The dataset is divided into training set (60% genes),
# validation set (15% genes), and test set (25% genes).

# limit dataset entries, useful to test models on a smaller dataset.
# set to -1 to disable
LIMIT_GENE_COUNT = 3000 * 400  # multiply by bin size

def load_rna_counts():
    esc_rna = pd.read_csv("./data/" + "RNA_counts_cleaned.txt", sep='\t')
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

def split_dataset():
    files = [f for f in listdir(DATA_DIR) if isfile(join(DATA_DIR, f)) and f.endswith(".all.csv")]

    rna_counts = load_rna_counts()

    for file in files:
        sample_name = file[0:len(file) - len(".all.csv")]
        data = pd.read_csv(join(DATA_DIR, file), header=None, names=["gene", "bin", "meth", "chr_acc"], sep='\t')
        total_cnt = data.shape[0]
        if LIMIT_GENE_COUNT != -1:
            total_cnt = min(data.shape[0], LIMIT_GENE_COUNT)
        train_loc = int(0.6 * total_cnt)
        while not data.iloc[train_loc]["bin"] == 0:
            train_loc += 1
        val_loc = int(0.75 * total_cnt)
        while not data.iloc[val_loc]["bin"] == 0:
            val_loc += 1

        train_set = data.iloc[:train_loc]
        val_set = data.iloc[train_loc:val_loc]
        test_set = data.iloc[val_loc:total_cnt]

        print(train_set.shape, val_set.shape, test_set.shape)

        train_set.to_csv(DATA_DIR + sample_name + ".train.csv", sep=",", header=False, index=False)
        val_set.to_csv(DATA_DIR + sample_name + ".valid.csv", sep=",", header=False, index=False)
        test_set.to_csv(DATA_DIR + sample_name + ".test.csv", sep=",", header=False, index=False)


split_dataset()
