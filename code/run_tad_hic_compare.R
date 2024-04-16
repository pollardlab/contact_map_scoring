# Run TADCompare and HiCcompare on genomic windows of ESC and HFF MicroC or HiC cool data

# Data needed: Two text files derived from cool files to compare at wanted resolution

library('TADCompare')
library('HiCcompare')


#### Fill out the following section ####

# Enter chromosomes to analyze
chroms = c('chr21', 'chr22')

# Enter resolutions to analyze
resolutions = c('2048', '4096', '8192', '10240', '20480')

# Enter prefix for input files. Input files should be Input file format: {input_file1}_{experiment}_{chrom}_{resolution}.txt
input_file1 = 'cool1'
input_file2 = 'cool2'

# Enter output directory. A new directory will be created within that separately for HiC compare and TADcompare results
output_dir = 'scores/'

# Enter prefix for output files. Output file names will be: {output_file}_{experiment}_{chrom}_{resolution}.txt
output_file = 'results_all_res' 


# Enter TRUE or FALSE if you want to or not run each method
run_HiCcompare_bool = TRUE
run_TADcompare_bool = TRUE





#### The following section will run HiCcompare and TADcompare for each experiment, chromosomes, and resolution provided ####


# HiC compare
run_HiCcompare <- function(cool_file_1, cool_file_2, output_file){
  
  hic_table <- create.hic.table(cool_file_1, cool_file_2)
  hic_table_norm = HiCcompare::hic_loess(hic_table)
  hic_table_result <- HiCcompare::hic_compare(hic_table_norm, Plot = F)
  write.table(hic_table_result, file = paste(output_dir, 'HiCcompare/', output_file, '.txt', sep = ''), 
              sep = '\t', col.names = T, quote = F, row.names = F)

}


# TADCompare
run_TADcompare <- function(cool_file_1, cool_file_2, output_file){

  cool_file_1_sparse <- HiCcompare::cooler2sparse(cool_file_1)
  cool_file_2_sparse <- HiCcompare::cooler2sparse(cool_file_2)
  diff_tads = TADCompare(cool_file_1_sparse, cool_file_2_sparse)
  write.table(diff_tads$TAD_Frame, file = paste(output_dir, 'TADcompare/', output_file, '.txt', sep = ''), 
              sep = '\t', quote = F, row.names = F)
}


# Run
set.seed(1)
for (res in resolutions){
  for (chrom in chroms){
      for (experiment in experiments){

        cool_file_1 <- read.table(paste(input_file1, '_', experiment, '_', chrom, '_', res,'.txt', sep = ''), 
                              sep = '\t', col.names = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'IF1'))
        cool_file_2 <- read.table(paste(input_file2, '_', experiment, '_', chrom, '_', res,'.txt', sep = ''), 
                              sep = '\t', col.names = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'IF1'))

        if run_HiCcompare_bool{
            run_HiCcompare(cool_file_1, cool_file_2, paste(output_file, experiment, chrom, res, sep = '_'))
            }
          
        if run_TADcompare_book{
            run_TADcompare(cool_file_1, cool_file_2, paste(output_file, experiment, chrom, res, sep = '_'))
            }

    }
  }
}





