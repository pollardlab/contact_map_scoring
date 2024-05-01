#!/usr/bin/env python3

import numpy as np
import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--bedfile", "-b", type=str, required=True)
parser.add_argument("--tadfirst", "-t", type=str, required=True)
parser.add_argument("--tadsecond", "-s", type=str, required=True)
parser.add_argument("--resolution", "-r", type=int, required=True)
parser.add_argument("--distance", "-d", type=int, required=True)
parser.add_argument( "--output","-o",type=str)

args = parser.parse_args()
bed_file = args.bedfile
tad_file1 = args.tadfirst
tad_file2 = args.tadsecond
res = args.resolution
ther = args.distance
outfile = args.output

def get_chrom_tad(tadfile):
    tads_dict = defaultdict(list)
    with open(tadfile,"r") as f:
        for line in f:
            if 'x1'  in line:
                continue
            info = line.strip().split('\t')
            chrom = 'chr'+info[0]
            if chrom in tads_dict:
                tads_dict[chrom].append((chrom,int(info[1]),int(info[2])))
            else:
                tads_dict[chrom] = [(chrom,int(info[1]),int(info[2]))]
    return tads_dict

def merge_adjacent_tads(tad_list,ther,res):
    merged_TAD = []
    tad_sorted = sorted(list(set(tad_list)))
    diff_index = np.flatnonzero(np.r_[True,np.diff(tad_sorted)>ther])
    if len(diff_index) == len(tad_sorted) or len(tad_sorted)<=1:
        for ts in tad_sorted:
            merged_TAD.append((ts,ts+res))
    else:        
        for i in range(len(diff_index)-1):
            if(diff_index[i+1]==(diff_index[i]+1)):
                merged_TAD.append((tad_sorted[diff_index[i]],tad_sorted[diff_index[i]]+res))
            else:
                merged_TAD.append((tad_sorted[diff_index[i]],tad_sorted[diff_index[i+1]-1]+res))
        
        if(diff_index[-1]==(len(tad_sorted)-1)):
            merged_TAD.append((tad_sorted[-1],tad_sorted[-1]+res))
        else:
            merged_TAD.append((tad_sorted[diff_index[-1]],tad_sorted[-1]+res))
    return merged_TAD
     
        
def get_tad_overlap_stats(bed_file,tad_s1,tad_s2,res,ther,outfile):
    fout = open(outfile, "w")
    fout.write('\t'.join(['chrom','start','end','H1hESC_TAD','HFFc6_TAD','H1hESC_overlap_TAD','HFFc6_overlap_TAD', \
                            'H1hESC_specific_TAD(loss)', 'HFFc6_specific_TAD(gain)','overlap_ratio','loss_ratio','gain_ratio','\n']))

    with open(bed_file,'r') as f:
        for line in f:
            reg = line.strip().split('\t')
            chrom=reg[0]
            start = int(float(reg[1]))
            end = int(float(reg[2]))
            tad_list1 = []
            tad_list2 = []
            for tr1 in tad_s1[chrom]:
                if ((tr1[1]>=start) and ((tr1[2]+res)<=end)):
                    tad_list1.append(tr1[1])
                    tad_list1.append(tr1[2])
            merged_tad_list1 = merge_adjacent_tads(tad_list1,ther,res)
            
            for tr2 in tad_s2[chrom]:
                if ((tr2[1]>=start) and ((tr2[2]+res)<=end)):
                    tad_list2.append(tr2[1])
                    tad_list2.append(tr2[2])
            merged_tad_list2 = merge_adjacent_tads(tad_list2,ther,res)
            
            overlap_tad_list1 = [x for x in merged_tad_list1 for y in merged_tad_list2 if (abs(y[0]-x[1])<=ther
                                 or abs(x[0]-y[1])<=ther or abs(x[0]-y[0])<=ther or abs(x[1]-y[1])<=ther)]
            overlap_tad_list2 = [y for x in merged_tad_list1 for y in merged_tad_list2 if (abs(y[0]-x[1])<=ther 
                                or abs(x[0]-y[1])<=ther or abs(x[0]-y[0])<=ther or abs(x[1]-y[1])<=ther)]
            tad_loss = list(set(merged_tad_list1) - set(overlap_tad_list1))
            tad_gain = list(set(merged_tad_list2) - set(overlap_tad_list2))
            merged_tad_1 = ';'.join([str(x) for x in merged_tad_list1]) 
            merged_tad_2 = ';'.join([str(x) for x in merged_tad_list2])
            overlap_tad_1 = ';'.join([str(x) for x in list(set(overlap_tad_list1))])
            overlap_tad_2 = ';'.join([str(x) for x in list(set(overlap_tad_list2))]) 
            tad_gain_str = ';'.join([str(x) for x in tad_gain])
            tad_loss_str = ';'.join([str(x) for x in tad_loss])
            max_tad = max(len(merged_tad_list1),len(merged_tad_list2))    
            overlap_ratio1 = len(set(overlap_tad_list1))/max_tad if max_tad > 0 else 1
            overlap_ratio2 = len(set(overlap_tad_list2))/max_tad if max_tad > 0 else 1
            overlap_ratio = overlap_ratio1 if overlap_ratio1>=overlap_ratio2 else overlap_ratio2 
            gain_ratio = len(tad_gain)/max_tad if max_tad > 0 else 0
            loss_ratio = len(tad_loss)/max_tad if max_tad > 0 else 0
            fout.write('\t'.join([chrom,reg[1],reg[2],merged_tad_1,merged_tad_2,overlap_tad_1,overlap_tad_2, tad_loss_str, \
                                tad_gain_str, str(overlap_ratio),str(loss_ratio),str(gain_ratio),'\n']))
                                #tad_gain_str, str(overlap_ratio1),str(overlap_ratio2),str(loss_ratio),str(gain_ratio),'\n']))


if __name__ == '__main__':
    tad_s1 = get_chrom_tad(tad_file1)
    tad_s2 = get_chrom_tad(tad_file2)
    get_tad_overlap_stats(bed_file,tad_s1,tad_s2,res,ther,outfile)
