import math

import pandas as pd
import os
import numpy as np

DATA_DIR = "./data/"
PROCESSED_DATA_DIR = "./processed_data/"
from itertools import combinations
import scipy.stats

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



def analyze_cell_data(cell_name):
    esc_rna_counts = load_rna_counts()
    data = pd.read_csv(PROCESSED_DATA_DIR + cell_name + ".all.csv", header=None, names=["gene", "bin", "meth", "chr_acc"], sep='\t')

    used_genes = ["ENSMUSG00000022652","ENSMUSG00000028164","ENSMUSG00000004665","ENSMUSG00000060499","ENSMUSG00000028184","ENSMUSG00000002076","ENSMUSG00000000365","ENSMUSG00000026970","ENSMUSG00000038733","ENSMUSG00000056515","ENSMUSG00000060499","ENSMUSG00000028184","ENSMUSG00000028184","ENSMUSG00000030672","ENSMUSG00000039329","ENSMUSG00000051176","ENSMUSG00000022652","ENSMUSG00000028164","ENSMUSG00000022652","ENSMUSG00000024968","ENSMUSG00000028164","ENSMUSG00000022221","ENSMUSG00000025480","ENSMUSG00000002489","ENSMUSG00000032525","ENSMUSG00000054889","ENSMUSG00000022652","ENSMUSG00000056515","ENSMUSG00000028164","ENSMUSG00000004665","ENSMUSG00000004665","ENSMUSG00000008393","ENSMUSG00000039395","ENSMUSG00000024030","ENSMUSG00000006134","ENSMUSG00000038742","ENSMUSG00000039329","ENSMUSG00000055148","ENSMUSG00000060461","ENSMUSG00000021255","ENSMUSG00000021255","ENSMUSG00000022652","ENSMUSG00000021255","ENSMUSG00000039329","ENSMUSG00000078952","ENSMUSG00000021255","ENSMUSG00000039329","ENSMUSG00000021255","ENSMUSG00000039329","ENSMUSG00000060499","ENSMUSG00000060461","ENSMUSG00000039209","ENSMUSG00000039329","ENSMUSG00000060499","ENSMUSG00000039209","ENSMUSG00000060461","ENSMUSG00000045699","ENSMUSG00000084762","ENSMUSG00000046761","ENSMUSG00000030491","ENSMUSG00000031758","ENSMUSG00000029106","ENSMUSG00000004665","ENSMUSG00000021255"]
    genes = data[data["gene"].isin(used_genes)]
    genes = genes.groupby(["gene"])

    # correlation among genes
    met_rates = []
    chr_rates = []
    expr_rates = []
    for name, group in genes:
        met_rate = group["meth"].sum()
        chr_rate = group["chr_acc"].sum()

        met_rate = (met_rate+1.0)
        chr_rate = (chr_rate+1.0)

        met_rates.append(met_rate)
        chr_rates.append(chr_rate)

        expr_rate = esc_rna_counts[name][cell_name[len("Cell_"):]]
        expr_rate = math.log10(expr_rate+1.0)
        expr_rates.append(expr_rate)

        print(name, met_rate, chr_rate, expr_rate)

    r, p = scipy.stats.pearsonr(met_rates, expr_rates)
    print(r ,p)
    r, p = scipy.stats.pearsonr(chr_rates, expr_rates)
    print(r ,p)



def analyze_data_across_celss(cell_names):
    esc_rna_counts = load_rna_counts()
    gene_dict = {}

    for cell_name in cell_names:

        data = pd.read_csv(PROCESSED_DATA_DIR + "Cell_" + cell_name + ".all.csv", header=None, names=["gene", "bin", "meth", "chr_acc"], sep='\t')

        genes = data.groupby(["gene"])

        for name, group in genes:
            met_rate = group["meth"].sum()
            chr_rate = group["chr_acc"].sum()

            met_rate = (met_rate+1.0)
            chr_rate = (chr_rate+1.0)

            expr_rate = esc_rna_counts[name][cell_name]
            expr_rate = math.log10(expr_rate+1.0)

            if name not in gene_dict:
                gene_dict[name] = []

            gene_dict[name].append((met_rate, chr_rate, expr_rate))

            # print(name, met_rate, chr_rate, expr_rate)

    gene_list = []
    for gene in gene_dict:
        met_rates = [ v[0] for v in gene_dict[gene] ]
        chr_rates = [ v[1] for v in gene_dict[gene] ]
        expr_rates = [v[2] for v in gene_dict[gene] ]
        r1, p1 = scipy.stats.pearsonr(met_rates, expr_rates)
        r2, p2 = scipy.stats.pearsonr(chr_rates, expr_rates)
        # if  r1 < -0.4 and p1 < 1.0:
        if r1 > 0.15:
            print(gene, r1, p1, r2, p2)
            gene_list.append(gene)

    print(gene_list)
    print(len(gene_list))


"""
for cell in ["Cell_ESC_H05", "Cell_ESC_B07", "Cell_ESC_H06", "Cell_ESC_H03"]:
    print(cell)
    analyze_cell_data(cell)
"""



cells = ["A02","A03","A04","A05","A06","A07","A09","B01","B02","B03","B04","B05","B06","B07","B09","C01","C02","C03","C04","C05","C07","C09","D01","D05","D06","D07","D09","E01","E02","E03","E05","E07","E09","F02","F03","F04","F05","F06","F07","F09","G01","G02","G03","G05","G06","G07","G09","H02","H03","H09"]
cells = ["ESC_" + v for v in cells]

analyze_data_across_celss(cells)

