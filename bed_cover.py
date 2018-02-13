#!/usr/bin/python3.6
from __future__ import print_function
import sys
import math
import argparse
import ntpath
import os

parser = argparse.ArgumentParser(description='Generate coverage info',
                                 epilog  = "usage example: bed_cov -cl 4 -cv 50 *_stat.bed")
parser.add_argument('files', metavar='files', type=str, nargs='+',
                    help='bed-files (example: *_stat.bed)')
parser.add_argument('-j', "--jaccard", dest="is_jaccard", action='store_true', default=False,
                    help='output jacard index instead of bed-file with ranges')
parser.add_argument('-cl', "--column", dest="column", metavar='column', type=int,
                    help='column to use (default: 4)', default=4)
parser.add_argument('-cv', "--coverage", dest="coverage", metavar='coverage (default: 50)', type=float,
                    help='coverage threshold', default=50)


args = parser.parse_args()
#print(args)

import glob

merged = {}

is_jaccard = args.is_jaccard
column = args.column - 1
coverage = args.coverage

if(len(args.files)==0):
    print("Files not found")
    exit(1)

files = []
[files.extend(glob.glob(file)) for file in args.files]

total_files = len(files)

for filename in files:
    with open(filename) as f:
        for line in f:
            temp = line.split('\t')
            temp[-1] = temp[-1].rstrip()
            curr_key = temp[0] + "\t" + temp[1] + "\t" + temp[2]
            if (curr_key not in merged):
                merged[curr_key] = []
            merged[curr_key].extend([float(v) for v in temp[3:]])

columns_total = len(next(iter(merged.items()))[1])
columns_in_file = columns_total / total_files

dNTP_all_genomic = 0
dNTP_all = 0
dNTP_used = 0

chr_old = 0
end_old = 0

if (is_jaccard):
    print("Total","Anchored","Covered", sep='\t')
else:
    print("chr", "start", "end", '\t'.join([ntpath.basename(file) for file in files]), sep='\t')
current_range = range(int(column), int(columns_total), int(columns_in_file))
for item in iter(merged.items()):

    if (is_jaccard):
        temp = item[0].split('\t')
        chr = temp[0]
        start = int(temp[1])
        end = int(temp[2])

        dNTP_all += (end - start)
        if(chr==chr_old):
            dNTP_all_genomic += ((end - start) + (start - end_old))
        else:
            dNTP_all_genomic += (end - start)

    if all(item[1][i] >= coverage for i in current_range):
        if(is_jaccard):
            dNTP_used += (end-start)
        else:
            print(item[0], '\t'.join( str(item[1][i]) for i in current_range), sep='\t')

    if (is_jaccard):
        chr_old = chr
        end_old = end
print(dNTP_all_genomic, dNTP_all, dNTP_used)
if (is_jaccard):
    print("A%T", "C%T", "C%A", sep='\t')
    print(round(float(dNTP_all)/float(dNTP_all_genomic),3),
          round(float(dNTP_used)/float(dNTP_all_genomic),3),
          round(float(dNTP_used)/float(dNTP_all),3), sep='\t')
