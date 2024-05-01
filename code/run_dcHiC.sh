#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
    Usage: run_dcHiC.sh [-C <cool_file1>] [-c <cool_file2>] [-r <region_file>] [-P <prefix1>]
                [-p <prefix2>] [-g <genome_size_file>] [-b <bin_size>] [-t <threads>] [-o <outdir>] [-d dchicdir]
    This script aims to run dcHiC from genome contact maps (.cool file) of two conditions.
    -C : The genome-wide contact map in .cool format from condition 1.
    -c : The genome-wide contact map in .cool format from condition 2.
    -r : The regions for map comparison.
    -P : The prefix of condition1.
    -p : The prefix of condition2.
    -g : The file of genome size.
    -b : bin size for contact map. Default is 2048.
    -t : The CPU threads for parallelized processing. Default is 1.
    -o : The output directory.
    -d : The directory of dcHiC.
    -h : Show usage help
EOF
    exit 1
}

while getopts :C:c:r:P:p:g:b:t:o:d:h opt; do
    case $opt in
        C) cool_file1=${OPTARG};;
        c) cool_file2=${OPTARG};;
        r) region_file=${OPTARG};;
        P) prefix1=${OPTARG};;
        p) prefix2=${OPTARG};;
        g) genome_size_file=${OPTARG};;
        b) bin_size=${OPTARG};;
        t) threads=${OPTARG};;
        o) outdir=${OPTARG};;
        d) dchicdir=${OPTARG};;
        h) usage;;
    esac
done

[ -z "$prefix1" ] && echo "Error! Please provide the prefix of condition1" && usage
[ -z "$prefix2" ] && echo "Error! Please provide the prefix of condition2" && usage
[ -z "$outdir" ] && echo "Error! Please provide the output directory" && usage
[ -z "$dchicdir" ] && echo "Error! Please provide the script directory of dcHiC" && usage
[ -z "$threads" ] && threads=1

if ! [[ -e "$cool_file1" && -e "$cool_file2" && -e "$region_file" ]]; then
      echo "Error! Two genome-wide contact maps, a region file and a dcHiC input file are needed." && usage 
else
    python ${dchicdir}/utility/preprocess.py -input cool -file $cool_file1 -genomeFile $genome_size_file -res $bin_size -prefix $prefix1
    python ${dchicdir}/utility/preprocess.py -input cool -file $cool_file2 -genomeFile $genome_size_file -res $bin_size -prefix $prefix2
    input_file=input.${prefix1}_${prefix2}_${bin_size}.dcHiC.txt
    region_file1=${prefix1}_${bin_size}_abs.bed
    region_file2=${prefix2}_${bin_size}_abs.bed
    cat > ${input_file} <<EOF
${prefix1}_${bin_size}.matrix	$region_file1	${prefix1}_${bin_size}	${prefix1}
${prefix2}_${bin_size}.matrix	$region_file2	${prefix2}_${bin_size}	${prefix2}
EOF
    Rscript ${dchicdir}/dchicf.r --file $input_file --pcatype cis --dirovwt T --cthread $threads --pthread $threads
    Rscript ${dchicdir}/dchicf.r --file $input_file --pcatype select --dirovwt T --genome hg38
    Rscript ${dchicdir}/dchicf.r --file $input_file --pcatype analyze --dirovwt T --diffdir $outdir
    python get_dcHiC_result.py --region_file $region_file --genome_size_file $genome_size_file --prefix1 $prefix1 --prefix2 $prefix2 --res $bin_size --datadir DifferentialResult/$outdir/pcQnm --outfile ${prefix1}_${prefix2}_${bin_size}_dcHiC_results.tsv
fi
