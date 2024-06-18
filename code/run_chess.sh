#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
    Usage: run_chess.sh [-C <cool_file1>] [-c <cool_file2>] [-R <region_comp_file>] 
                [-b <bin_size>] [-t <threads>] [-o <output_prefix>]
    This script aims to run chess from genome contact maps (.cool file) or predicted contact maps (.h5 file) of two conditions. Either genome contact maps or predicted contact maps are needed.
    -C : The genome-wide contact map in .cool format from condition 1.
    -c : The genome-wide contact map in .cool format from condition 2.
    -R : The regions for map comparison in bedpe format: chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2.
    -b : bin size for contact map. Default is 2048.
    -t : The CPU threads for parallelized processing. Default is 8.
    -o : The prefix of output files.
    -d : The output file(s) directory.
    -h : Show usage help
EOF
    exit 1
}

while getopts :C:c:R:b:t:o:h opt; do
    case $opt in
        C) cool_file1=${OPTARG};;
        c) cool_file2=${OPTARG};;
        R) region_comp_file=${OPTARG};;
        b) bin_size=${OPTARG};;
        t) threads=${OPTARG};;
        o) output_prefix=${OPTARG};;
        d) output_dir=${OPTARG};;
        h) usage;;
    esac
done

[ -z "$output_prefix" ] && echo "Error! Please provide the prefix of output files" && usage
[ -z "$bin_size" ] && bin_size=2048
[ -z "$threads" ] && threads=8
if ! [[ -e "$cool_file1" && -e "$cool_file2" && -e "$region_comp_file" ]]; then
    echo "Error! Two genome-wide contact maps and a region file in bedpe format are needed for experimental data." && usage 
else
    if [[ $cool_file1 == *.mcool ]]; then
        chess sim -p $threads ${cool_file1}::resolutions/${bin_size} ${cool_file2}::resolutions/${bin_size} $region_comp_file ${output_file}_details.tsv
    else
        chess sim -p $threads $cool_file1 $cool_file2 $region_comp_file ${output_file}_details.tsv
    fi
    cat ${output_file}.tsv |grep -v 'ID'| paste $region_comp_file - |tr " " "\t" |cut -f1,2,3,12,13|awk 'BEGIN{print "chrom\tstart\tend\tSN\tssim"}; {print $0}' > ${output_dir}${output_file}.tsv  
fi
