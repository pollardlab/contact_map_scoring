#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--bedfile", "-b", type=str, required=True)
parser.add_argument("--diffloopfile", "-d", type=str, required=True)
parser.add_argument( "--output","-o",type=str)

args = parser.parse_args()
bed_file = args.bedfile
diff_loop_file = args.diffloopfile
outfile = args.output

def get_overlap_count(bed_file,diff_loop_file,outfile):
    diff_loops = open(diff_loop_file,"r").readlines()
    fout = open(outfile, "w")
    with open(bed_file,'r') as f:
        for line in f:
            reg = line.strip().split('\t')
            chrom=reg[0]
            start = int(float(reg[1]))
            end = int(float(reg[2]))
            overlap_count =0
            for dl in diff_loops:
                if('CHR1' in dl):
                    continue
                loop_chrom = dl.strip().split('\t')[0]
                loop_start = dl.strip().split('\t')[1]
                loop_end = dl.strip().split('\t')[5]
                if((chrom==loop_chrom) and (int(loop_start)>start) and (int(loop_end)<end)):
                    overlap_count+=1
            fout.write(chrom+"\t"+reg[1]+"\t"+reg[2]+"\t"+str(overlap_count)+"\n")


if __name__ == '__main__':
    get_overlap_count(bed_file,diff_loop_file,outfile)
