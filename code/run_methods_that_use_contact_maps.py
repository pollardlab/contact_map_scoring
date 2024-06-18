#!/usr/bin/env python

import pandas as pd
import cooler
import cooltools
from cooltools.lib.numutils import observed_over_expected
from cooltools.lib.numutils import adaptive_coarsegrain
from cooltools.lib.numutils import interpolate_bad_singletons
from cooltools.lib.numutils import interp_nan, set_diag
from cooltools.lib.plotting import *


########################################
#### Fill out the following section ####


# Enter file name for tab-delimited text file with genomic windows to score, that includes columns: chrom, start, end
file_name = '../data/experimental_maps/example_input/example_DEG_windows.bed'


# Enter file names for each experiment type. Comment out experiments you are not including.
experiment_files = {
    
    'MicroC':['../data/experimental_maps/ESC_MicroC.mcool', 
              '../data/experimental_maps/HFF_MicroC.mcool']
#     ,
#     'HiC':['cool_file1.mcool',
#            'cool_file2.mcool']
    
                   }


# Enter window sizes included and their resolution. Window size should be equal to end-start
window_sizes = {
#     '100kb':[100000, 2048]
#     ,
                
    '1Mb':[1000000, 2048]
#     ,
                
#     '10Mb':[10000000, 20480]
}


# Enter output directory.
output_dir = '../data/experimental_maps/example_output/'


# Enter output file name
output_file = 'methods_that_use_contact_maps_scores' 


# Choose scores to apply by commenting out ones to skip
score_types = [
    'correlation', # Spearman correlation 
       'mse', # Mean squared error
       'ssim_map', # Structural similarity index measure
       'contact_directionality_corr', # Contact directionality (corr)
       'distance_enrichment_corr', # Distance enrichment (corr)
       'eigenvector_corr', # Eigenvector (corr)
       'insulation_corr', # Insulation (corr)
       'insulation_mse', # Insulation (mse)
       'loops', # Loops
       'scc', # Stratum-adjusted correlation coefficient
       'tads', # TADs
       'triangle_corr', # Triangle (corr)
       'triangle_mse' # Triangle (mse)
              ]



#####################################################################################################
#### The following section will run the selected methods for each experiment and window provided ####


# Create output data frame

windows = pd.read_csv(file_name, sep = '\t')

chroms = np.unique(windows.chrom)

windows['window_id'] = windows[['chrom', 'start', 'end']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
windows['window'] = [x-y for x,y in zip(windows.end, windows.start)]

windows = windows[['window_id', 'chrom', 'start', 'end', 'window']]

# Duplicate dataset for each experiment
scores = pd.DataFrame()
for experiment in experiment_files.keys():
    scores_ = windows.copy()
    scores_['experiment'] = experiment
    scores = pd.concat([scores, scores_], axis = 0)




def get_experimental_map(myseq_str, genome_hic_cool):

    num_counts= np.sum(genome_hic_cool.matrix(balance=False).fetch(myseq_str))
    seq_hic_obs = genome_hic_cool.matrix(balance=True).fetch(myseq_str)
    
    seq_hic_smoothed =  adaptive_coarsegrain(
                     seq_hic_obs,  
                     genome_hic_cool.matrix(balance=False).fetch(myseq_str),  
                     cutoff=3, max_levels=8)
    seq_hic_nan = np.isnan(seq_hic_smoothed)
    seq_hic_obsexp = observed_over_expected(seq_hic_smoothed, ~seq_hic_nan)[0]
    seq_hic_obsexp = np.log(seq_hic_obsexp)
    seq_hic_obsexp = np.clip(seq_hic_obsexp,-2,2)
    seq_hic_obsexp_init = np.copy(seq_hic_obsexp)
    seq_hic_obsexp = interp_nan(seq_hic_obsexp)
    seq_hic_obsexp = np.nan_to_num(seq_hic_obsexp)
    seq_hic = np.clip(seq_hic_obsexp,-2,2)
    for i in [-1,0,1]: set_diag(seq_hic, 0,i)
        
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    kernel = Gaussian2DKernel(x_stddev=1,x_size=5)

    seq_hic = convolve(seq_hic, kernel)
    
    return seq_hic




# # # # # # # # # # # # # #
# Scoring functions 


import scoring


corr_scores = ['correlation',  
               'ssim_map', 
               'contact_directionality_corr', 
               'distance_enrichment_corr', 
               'eigenvector_corr', 
               'insulation_corr',  
               'loops', 
               'tads',
               'triangle_corr']

class scoring_map_methods:
    
    '''
    Define the function necessary for calculating disruption scores based on map-motivated scoring methods. 
    
    '''
    
    def __init__(self, map_a, map_b):
        self.map_a = map_a
        self.map_b = map_b

        
    # Matrix methods
    
    def correlation(self): # spearman correlation
        return scoring.correlation(self.map_a, self.map_b)   
    
    def mse(self): # mean squared error
        return scoring.mse(self.map_a, self.map_b)
    
    def ssim_map(self): # structural similarity index
        return scoring.ssim_map(self.map_a, self.map_b)
    
    
    # Contact map methods
    
    def contact_directionality_corr(self): # DI_track
        return scoring.vectorMethodToScalar(scoring.contact_directionality_track, self.map_a, self.map_b, finalCompMetric = 'corr')    
    
    def distance_enrichment_corr(self): # decay_track
        return scoring.vectorMethodToScalar(scoring.distance_enrichment_track, self.map_a, self.map_b, finalCompMetric = 'corr')   
    
    def eigenvector_corr(self): # contact_pca_track
        return scoring.vectorMethodToScalar(scoring.eigenvector_track, self.map_a, self.map_b, finalCompMetric = 'corr')

    def insulation_corr(self): # insulation_track
        return scoring.vectorMethodToScalar(scoring.insulation_track, self.map_a, self.map_b, finalCompMetric = 'corr')
      
    def insulation_mse(self): # insulation_track
        return scoring.vectorMethodToScalar(scoring.insulation_track, self.map_a, self.map_b, finalCompMetric = 'mse')
        
    def loops(self):
        return scoring.Loops(self.map_a, self.map_b)['overlap_ratio']   
    
    def scc(self): # stratum adjusted correlation
        return scoring.scc(self.map_a, self.map_b)
    
    def tads(self): # insulation_track
        return scoring.TAD(self.map_a, self.map_b)['overlap_ratio']
    
    def triangle_corr(self): # triangle_track
        return scoring.vectorMethodToScalar(scoring.triangle_track, self.map_a, self.map_b, finalCompMetric = 'corr')

    def triangle_mse(self): # triangle_track
        return scoring.vectorMethodToScalar(scoring.triangle_track, self.map_a, self.map_b, finalCompMetric = 'mse')
    
    

    



# # # # # # # # # # # # # #
# Run


# Get scores

for window_size in window_sizes.keys():
    
    for experiment in experiment_files.keys():
        
        for chrom in chroms:
        
            cool1 = cooler.Cooler(f'{experiment_files[experiment][0]}::/resolutions/{window_sizes[window_size][1]}')
            cool2 = cooler.Cooler(f'{experiment_files[experiment][1]}::/resolutions/{window_sizes[window_size][1]}')

            for i in scores[(scores.chrom == chrom) & 
                            (scores.experiment == experiment) &
                            (scores.window == window_sizes[window_size][0])].index:

                start = scores.loc[i,'start']
                end = scores.loc[i,'end']

                map1 = get_experimental_map(f'{chrom}:{start}-{end}', cool1)
                map2 = get_experimental_map(f'{chrom}:{start}-{end}', cool2)


                for score in score_types:
                    scores.loc[i,score] = getattr(scoring_map_methods(map1.copy(),map2.copy()), score)()
        

    
windows_to_score.to_csv(f'{output_dir}{output_file}', sep ='\t', index = False)
    
    
    
    
    




