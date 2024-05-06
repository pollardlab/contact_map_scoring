#!/bin/bash

########################################
#### Fill out the following section ####

cool_file1='../data/experimental_maps/ESC_MicroC.mcool'
cool_file2='../data/experimental_maps/HFF_MicroC.mcool'
window_file='../data/experimental_maps/example_DEG_windows_noheader.bed'
sample_prefixs=("ESC" "HFF")

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
	selfish -f1 $cool_file1 -f2 $cool_file2 -ch ${arrchr[$i]} -r ${arrres[$j]} -t 0.05 -o ./${sample_prefixs[0]}_${sample_prefixs[1]}_selfish_${arrchr[$i]}_${arrres[$j]}.tsv
       	window_file_prefix="${window_file%.*}"
	awk -F "\t" -v chrom="${arrchr[$i]}" '{if($1==chrom){print $0}}' $window_file > ${window_file_prefix}_${arrchr[$i]}.bed
        python get_overlap_count_selfish.py -b ${window_file_prefix}_${arrchr[$i]}.bed -d ${sample_prefixs[0]}_${sample_prefixs[1]}_selfish_${arrchr[$i]}_${arrres[$j]}.tsv -o ${window_file_prefix}_${arrchr[$i]}_${arrres[$j]}_selfish.tsv
    done
done
