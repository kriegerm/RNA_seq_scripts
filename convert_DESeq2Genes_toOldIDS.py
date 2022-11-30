#!/usr/bin/env python3

import sys
from Bio import SeqIO
import csv

#convert_DESeq2Genes_toOldIDS.py
#LAST UPDATE: 11/30/22
#It's very annoying how DESEq2 usually gives a mixture of locus tags and gene names, and it's even more annoying when you're trying to use the OLD locus tags for something (like kEGG analysis!)
#So this code takes in a csv of DEGs with the "Gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj" header.
#It'll take in a NewGenbank file with the fancy new identifiers that they have for species like S. mutans, and then give you a csv
#With all the important info from the DEG results as well as the new locus tag.

def DEGS_to_KEGG(DEGS_csv, NewGenbank, outputcsv):

    old_locus_tags_list = []
    log2FC_list = []
    pval_list = []
    padj_list = []

    with open(DEGS_csv, 'r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # skip the headers
        for line in reader:
            name = line[0]
            log2FC = line[2]
            pval = line[5]
            padj = line [6]

            old_locus_tags_list.append(name.strip())
            log2FC_list.append(log2FC.strip())
            pval_list.append(pval.strip())
            padj_list.append(padj.strip())


    #Iterate through the list of locus tags
    results_list = []
    for item in old_locus_tags_list:
        count = 0

        #Open the gb file and read it in using SeqIO
        with open(NewGenbank) as input_handle:
            for gb_record in SeqIO.parse(input_handle, "genbank"):

                #Look through each record in the gb
                for f in gb_record.features:
                    new_tag = str(f.qualifiers.get('locus_tag', [])) #This is the NEW locus tag
                    gene = str(f.qualifiers.get('gene', [])) #This is the gene name, if it's there
                    new_tag = new_tag.lstrip("['").rstrip("']'")
                    gene = gene.lstrip("['").rstrip("']'")

                    #And see if the locus tag of the entry matches the item
                    if item == new_tag:
                        #We only want a single entries per locus tag (becuase there are both gene and CDS entires)
                        if count == 0:

                            #So once we find one, we'll add 1 to the count and stop this loop so the script can move onto the next
                            count = count + 1
                            original_tag = list(f.qualifiers.get('old_locus_tag', [])) #This is the original (or old) locus tag

                            #Sometimes there are more than one original locus tag listed in the GB file (like SMU.01 and SMU_01), so depending on how many are there we only want to keep one result per gene
                            if len(original_tag) == 2:
                                results_list.append(original_tag[1])

                            if len(original_tag) == 1:
                                results_list.append(original_tag[0])

                            if len(original_tag) == 0:
                                results_list.append("NA")

                    #If we went through the entire gb and didn't find a hit for the locus tag, then the DEG results probably have the gene name (like "gyrA") on there instad of the tag
                    #We only want a single entries per locus tag (becuase there are both gene and CDS entires)
                    if count == 0:

                        #And see if the locus tag of the entry matches the item
                        if item == gene:

                            #So once we find one, we'll add 1 to the count and stop this loop so the script can move onto the next
                            count = count + 1
                            original_tag = list(f.qualifiers.get('old_locus_tag', []))
                            #results_list.append(original_tag)

                            #Sometimes there are more than one original locus tag listed in the GB file (like SMU.01 and SMU_01), so depending on how many are there we only want to keep one result per gene
                            if len(original_tag) == 2:
                                results_list.append(original_tag[1])

                            if len(original_tag) == 1:
                                results_list.append(original_tag[0])

                            if len(original_tag) == 0:
                                results_list.append("NO OLDER GENE NAME AVALIABLE")

    #Checking to make sure everything was stored properly in lists now that we are done looping
    if len(results_list) != len(pval_list) != len(padj_list) != len(old_locus_tags_list):
        print("WARNING: LENGTHS OF LISTS DO NOT MATCH!!!")

    #Write a header
    with open(outputcsv, 'a') as fh:
            writer = csv.writer(fh)
            header = ['Original Gene_Name', 'locus_tag', 'Log2FC', 'pvalue', 'padj']
            writer.writerow(header)

    #Write the data back to a csv file as  your results
    for i in range(len(results_list)):
        row = []
        row = [old_locus_tags_list[i], results_list[i], log2FC_list[i], pval_list[i], padj_list[i]]

        with open(outputcsv, 'a') as fh:
            writer = csv.writer(fh)
            writer.writerow(row)

#This little block just handles how many input variables you give it to make sure you have the correct number of arguments
if __name__ == '__main__':
    if len(sys.argv) == 4:
         DEGS_to_KEGG(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
         print("Usage: gb, locus tags csv, outputcsv")
         sys.exit(0)
