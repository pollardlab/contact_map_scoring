#!/usr/bin/env python3

import os
import sys, subprocess, time
import numpy as np
import pandas as pd
import argparse
import cooler

parser = argparse.ArgumentParser()
parser.add_argument("--regfile", "-i", type=str, required=True)
parser.add_argument("--coolfile1", "-f", type=str, required=True)
parser.add_argument("--coolfile2", "-s", type=str, required=True)
parser.add_argument("--sample_name1", "-a", type=str, required=True)
parser.add_argument("--sample_name2", "-b", type=str, required=True)
parser.add_argument("--resolution", "-r", type=int, required=True)
parser.add_argument("--chrom", "-c", type=str, required=True)
parser.add_argument("--window", "-w", type=int, required=True)
parser.add_argument("--processes", "-n", type=int, required=True)
parser.add_argument( "--outdir","-o",type=str)

args = parser.parse_args()
reg_file = args.regfile
cool_file1 = args.coolfile1
cool_file2 = args.coolfile2
sample_name1 = args.sample_name1
sample_name2 = args.sample_name2
res = args.resolution
chrom = args.chrom
window_size = args.window
num_processes = args.processes
outdir = args.outdir

## From basenji
def exec_par(cmds, max_proc=None, verbose=False):
    total = len(cmds)
    finished = 0
    running = 0
    p = []

    if max_proc == None:
        max_proc = len(cmds)

    if max_proc == 1:
        while finished < total:
            if verbose:
                print(cmds[finished], file=sys.stderr)
            op = subprocess.Popen(cmds[finished], shell=True)
            os.waitpid(op.pid, 0)
            finished += 1

    else:
        while finished + running < total:
            # launch jobs up to max
            while running < max_proc and finished+running < total:
                if verbose:
                    print(cmds[finished+running], file=sys.stderr)
                p.append(subprocess.Popen(cmds[finished+running], shell=True))
                #print 'Running %d' % p[running].pid
                running += 1

            # are any jobs finished
            new_p = []
            for i in range(len(p)):
                if p[i].poll() != None:
                    running -= 1
                    finished += 1
                else:
                    new_p.append(p[i])

            # if none finished, sleep
            if len(new_p) == len(p):
                time.sleep(1)
            p = new_p

        # wait for all to finish
        for i in range(len(p)):
            p[i].wait()

def run_HiC1Dmetrics_target_matrix(df,res,genome_hic_cool_1,genome_hic_cool_2,sample_name1,sample_name2,num_processes,outdir):
    hic1d_jobs = []

    num_regs = df.shape[0]
    for i in range(num_regs):
        chrom = df.iloc[i]['chrom']
        start = df.iloc[i]['window_start']
        end = df.iloc[i]['window_end']
        mseq_str = '%s:%d-%d' % (chrom, start, end)
        seq_hic_raw_s1 = genome_hic_cool_1.matrix(balance=True).fetch(mseq_str)
        seq_hic_raw_s1_row = np.r_[[(start//res)*res+res*np.array(range(seq_hic_raw_s1.shape[0]))],seq_hic_raw_s1]
        seq_hic_raw_s1_row_column = np.column_stack((np.append(np.nan,(start//res)*res+res*np.array(range(seq_hic_raw_s1.shape[0]))),seq_hic_raw_s1_row)) 
        
        seq_hic_raw_s2 = genome_hic_cool_2.matrix(balance=True).fetch(mseq_str)
        seq_hic_raw_s2_row = np.r_[[(start//res)*res+res*np.array(range(seq_hic_raw_s2.shape[0]))],seq_hic_raw_s2]
        seq_hic_raw_s2_row_column = np.column_stack((np.append(np.nan,(start//res)*res+res*np.array(range(seq_hic_raw_s2.shape[0]))),seq_hic_raw_s2_row))
        
        try:
            os.makedirs(f'{outdir}/{chrom}/{res}')
        except:
            pass
        np.savetxt(f'{outdir}/{chrom}/{res}/{sample_name1}_{res}_{i}_regions.txt',seq_hic_raw_s1_row_column,delimiter='\t')
        np.savetxt(f'{outdir}/{chrom}/{res}/{sample_name2}_{res}_{i}_regions.txt',seq_hic_raw_s2_row_column,delimiter='\t')
        
        if(res==20480):
            param = 102400
            param_DLR = 307200
        else:
            param = 10240
            param_DLR = 30720
        for metric in ['ISC','CIC','SSC','deltaDLR','CD']:
            cmd = f'h1d two {metric}'
            cmd += f' {outdir}/{chrom}/{res}/{sample_name1}_{res}_{i}_regions.txt'
            cmd += f' {outdir}/{chrom}/{res}/{sample_name2}_{res}_{i}_regions.txt'
            if(metric=='deltaDLR'):
                cmd += f' {res} {chrom} --datatype matrix -p {param_DLR}'
            elif(metric=='CD'):
                cmd += f' {res} {chrom} --datatype matrix'
            else:
                cmd += f' {res} {chrom} --datatype matrix -p {param}'
            cmd += f' {sample_name1}_vs_{sample_name2}_{i}_{metric}'
            hic1d_jobs.append(cmd)

    exec_par(hic1d_jobs, num_processes, verbose=True)

def extract_HiC1Dmetrics_results(df,res,sample_name1,sample_name2,outdir):
    metrics = ['ISC','CIC','SSC','deltaDLR','CD']
    num_regs = df.shape[0]
    
    for metric in metrics:
        metric_value_list = []
        for i in range(num_regs):
            chrom = df.iloc[i]['chrom']
            data_df = pd.read_table(f'{outdir}/{chrom}/{res}/{sample_name1}_vs_{sample_name2}_{i}_{metric}.bedGraph',
                                           header=None,sep='\t',names=['chrom','start','end','value'])
            metric_value = data_df['value'].abs().mean()
            metric_value_list.append(metric_value)
        df[metric] = metric_value_list
    df.to_csv(f'{outdir}/{sample_name1}_vs_{sample_name2}_comp_HiC1Dmetrics_results_{res}.tsv',sep='\t',index=False)  
            
if __name__ == '__main__':
    regs = pd.read_table(reg_file,header=0,sep='\t')
    regs_sub = regs[(regs['window']==window_size) & (regs['chrom']==chrom)]
    run_HiC1Dmetrics_target_matrix(regs_sub,res,cool_file1,cool_file2,sample_name1,sample_name2,num_processes,outdir)
    extract_HiC1Dmetrics_results(regs_sub,res,sample_name1,sample_name2,outdir)

           

