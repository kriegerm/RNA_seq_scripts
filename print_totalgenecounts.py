#!/usr/bin/env python3

##I didn't fully finish this code, I just wrote it quickly. It helps to run the code in terminal
#as $ ./script.py val1 val2 val3 > output.txt so you can save all your terminal output.
#Just get your counts of all the genes you'd input into DESeq2, and also a list of gene names as a
#.txt file that you want ot print out expresison values for. This will just look through the counts csv
#and print out the gene counts for the genes in the list. Probably a good idea to double check that
#the same number of gene sare in both lists in the same order!

import csv
import sys

def printlines(counts_csv, listofgenes):
    list_of_genes = []
    with open(listofgenes, 'r') as fh:
        for line in fh:
            list_of_genes.append(line.strip('\n'))

    counts_dict = {}
    with open(counts_csv, 'r') as fh:
        fhcsv = csv.reader(fh, delimiter=',')
        field_names_list = next(fhcsv)
        for entry in fhcsv:
            templist = []
            for line in entry:
                line.lstrip("''")
                line.rstrip("''")
                templist.append(line)
            counts_dict[templist[0]] = templist[1:]

    for entry in list_of_genes:
        if entry in counts_dict:
            print(entry, counts_dict[entry])


if __name__ == '__main__':
    if len(sys.argv) == 3:
        printlines(sys.argv[1], sys.argv[2])
    else:
        print("CHECK")
        sys.exit(0)
