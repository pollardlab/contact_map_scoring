#!/bin/bash

########################################
#### Fill out the following section ####

hic_file1='../data/experimental_maps/ESC_MicroC.hic'
hic_file2='../data/experimental_maps/HFF_MicroC.hic'
window_file='../data/experimental_maps/example_input/example_DEG_windows_noheader.bed'
sample_prefixs=("ESC" "HFF")
output_file='../data/experimental_maps/example_output/arrowhead.tsv'

# Insert chromosomes to analyze here, with quotes and space-separated
declare -a arrchr=("chr22")

# Insert resolutions to analyze here, space-separated
declare -a arrres=(2048)



################################################
#### The following section will run Arrowhead and get the TAD overlap between two samples####

for res in ${arrres[@]}; do
	java -Xmx10g -jar juicer_tools.jar arrowhead -m 2000 -r ${res} -k KR hic_file1 ${sample_prefixs[0]}_${res}_arrowhead
	java -Xmx10g -jar juicer_tools.jar arrowhead -m 2000 -r ${res} -k KR hic_file2 ${sample_prefixs[1]}_${res}_arrowhead
	for chr in ${arrchr[@]}; do
		window_file_prefix="${window_file%.*}"
		awk -F "\t" -v chrom="$chr" '{if($1==chrom){print $0}}' $window_file > ${window_file_prefix}_${chr}.bed
		python get_overlap_regions_arrowhead.py -b ${window_file_prefix}_${chr}.bed -t ${sample_prefixs[0]}_${res}_arrowhead/${res}_blocks.bedpe -s ${sample_prefixs[1]}_${res}_arrowhead/${res}_blocks.bedpe -r $res -d 10240 -o ${output_file}
	done
done

