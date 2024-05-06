#!/usr/bin/env python3

import io
import sys
import numpy as np
import pandas as pd
import scoring
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--region_file", help='genomic coordinates of the regions for pairwise comparison', required=True)
    parser.add_argument("--genome_size_file", help='chromosize size of the genome', required=True)
    parser.add_argument("--prefix1", help="prefix of results for condition 1", required=True)
    parser.add_argument("--prefix2", help="prefix of results for condition 2", required=True)
    parser.add_argument("--crop_length", type=int, default=0, help='crop bp off from each end of the genomic regions')
    parser.add_argument("--res", type=int, default=2048, help='the resolution of the predicted contact maps')
    parser.add_argument( "--datadir", help="directory of the dcHiC results", required=True)
    parser.add_argument( "--outfile", help="mse and pearson correlation of dcHiC compartmentalization", required=True)
    args = parser.parse_args()
    
    fout = open(args.outfile,'w')
    fout.write('\t'.join(['chrom','start','end','dcHiC_chrom','dcHiC_start','dcHiC_end','dcHiC_mse','dcHiC_pearson','dcHiC_spearman','\n']))
    dcHiC_Qnm = {}
    for line in open(args.genome_size_file):
        chrsize = line.strip().split()
        try:
            dcHiC_Qnm[chrsize[0]] = pd.read_table(f'{args.datadir}/intra_sample_{chrsize[0]}_combined.pcQnm.bedGraph',header=0,sep='\t')
        except FileNotFoundError as e:
            print(f'File intra_sample_{chrsize[0]}_combined.pcQnm.bedGraph not found!')
            pass
    
    for index,line in enumerate(open(args.region_file)):
        seq = line.strip().split()
        chrom = seq[0]
        org_start = int(seq[1])
        org_end = int(seq[2])
        start = (int(seq[1])//args.res)*args.res+args.crop_length
        end = (int(seq[2])//args.res)*args.res-args.crop_length
        try:
            Qnm_df = dcHiC_Qnm[chrom]
            pc_arr_1 = get_pc_array(chrom,start,end,args.res,f'{args.prefix1}_{str(args.res)}',Qnm_df)
            pc_arr_2 = get_pc_array(chrom,start,end,args.res,f'{args.prefix2}_{str(args.res)}',Qnm_df)
            print(pc_arr_1)
            #print(pc_arr_2)
            spearman = scoring.spearman_1D(pc_arr_1, pc_arr_2)
            pearson = scoring.pearson_1D(pc_arr_1, pc_arr_2)
            mse = scoring.mse_1D(pc_arr_1, pc_arr_2)
            fout.write('\t'.join([str(s) for s in [chrom,org_start,org_end,chrom,start,end,mse,pearson,spearman,'\n']]))
        except KeyError:
            print(f'dcHiC result for {chrom} does not exist.')

def get_pc_array(chrom,start,end,res,col,df):
    pc_arr = []
    for coor in range(start,end,res):
        try:
            rloc = df.index[(df['chr']==chrom) & (df['start']==coor)][0]
            rv = df.iloc[rloc][col]
            pc_arr.append(rv)
        except:
            pc_arr.append(0)
    return np.array(pc_arr)

if __name__ == '__main__':
    main() 
