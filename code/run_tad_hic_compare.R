# Run TADCompare and HiCcompare on genomic windows of ESC and HFF MicroC or HiC cool data

# Data needed: Two text files derived from cool files to compare at wanted resolution

library('TADCompare')
library('HiCcompare')


#### Fill out the following section ####

# Enter chromosome to run on.
chrom = 'chr21'

# Enter input files.
input_file1 = '../data/experimental_maps/example_input/ESC_MicroC_2048_chr21.txt'
input_file2 = '../data/experimental_maps/example_input/HFF_MicroC_2048_chr21.txt'

# Enter output directory. 
output_dir = '../data/experimental_maps/example_output/'


# Enter TRUE or FALSE if you want to or not run each method
run_HiCcompare_bool = TRUE
run_TADcompare_bool = TRUE



#### The following section will run HiCcompare and TADcompare for each experiment, chromosomes, and resolution provided ####


# HiC compare
run_HiCcompare <- function(cool_file_1, cool_file_2, output_file){
  
  hic_table <- create.hic.table(cool_file_1, cool_file_2)
  hic_table_norm = HiCcompare::hic_loess(hic_table)
  hic_table_result <- HiCcompare::hic_compare(hic_table_norm, Plot = F)
  write.table(hic_table_result, file = paste(output_dir, output_file, sep = ''), 
              sep = '\t', col.names = T, quote = F, row.names = F)

}


# TADCompare
run_TADcompare <- function(cool_file_1, cool_file_2, output_file){

  cool_file_1_sparse <- HiCcompare::cooler2sparse(cool_file_1)
  cool_file_2_sparse <- HiCcompare::cooler2sparse(cool_file_2)
  diff_tads = TADCompare(cool_file_1_sparse, cool_file_2_sparse)
  write.table(diff_tads$TAD_Frame, file = paste(output_dir, output_file, sep = ''), 
              sep = '\t', quote = F, row.names = F)
}


# Run
set.seed(1)

cool_file_1 <- read.table(input_file1, 
                      sep = '\t', col.names = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'IF1'))
cool_file_2 <- read.table(input_file2, 
                      sep = '\t', col.names = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'IF1'))

if run_HiCcompare_bool{
    run_HiCcompare(cool_file_1, cool_file_2, paste('HiCcompare_', chrom, '.txt', sep = ''))
    }
  
if run_TADcompare_book{
    run_TADcompare(cool_file_1, cool_file_2, paste('TADcompare_', chrom, '.txt', sep = ''))
    }





