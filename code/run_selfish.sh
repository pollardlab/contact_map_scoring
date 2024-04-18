#!/bin/bash

########################################
#### Fill out the following section ####

cool_file1 = '../data/experimental_maps/ESC_MicroC.mcool'
cool_file2 = '../data/experimental_maps/HFF_MicroC.mcool'

# Insert chromosomes to analyze here, with quotes and space-separated
declare -a arrchr=("chr22")

# Insert resolutions to analyze here, space-separated
declare -a arrres=(2048)



################################################
#### The following section will run Selfish ####


chrlength=${#arrchr[@]}
reslength=${#arrres[@]}
for (( i=0; i<${chrlength}; i++ ));
do
    for(( j=0; j<${reslength}; j++ ));
    do
	selfish -f1 $cool_file1 -f2 $cool_file2 -ch ${arrchr[$i]} -r ${arrres[$j]} -t 0.05 -o ./H1ESC_HFFc6_selfish_${arrchr[$i]}_${arrres[$j]}.tsv
        if [ -f H1_HFF_DEG_gene_regions_${arrchr[$i]}.bed ]; then
            python get_overlap_count_modify.py -b H1_HFF_DEG_gene_regions_${arrchr[$i]}.bed -d H1ESC_HFFc6_selfish_${arrchr[$i]}_${arrres[$j]}.tsv -o H1_HFF_DEG_gene_regions_${arrchr[$i]}_${arrres[$j]}_selfish.tsv
        else
            cat /pollard/data/projects/skuang/Scoring/H1_HFF_DEG/chess_MicroC/H1_HFF_DEG_gene_regions_comparison.bedpe |awk -v chrom="${arrchr[$i]}" -F "\t" '{if($1==chrom){print $1"\t"$2"\t"$3}}' > H1_HFF_DEG_gene_regions_${arrchr[$i]}.bed
            python get_overlap_count_modify.py -b H1_HFF_DEG_gene_regions_${arrchr[$i]}.bed -d H1ESC_HFFc6_selfish_${arrchr[$i]}_${arrres[$j]}.tsv -o H1_HFF_DEG_gene_regions_${arrchr[$i]}_${arrres[$j]}_selfish.tsv
        fi
    done
done
