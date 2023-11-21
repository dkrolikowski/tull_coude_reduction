""" Functions to measure radial velocities from extracted echelle spectra.

Created by DMK on 11/20/2023
Last updated by DMK on 11/20/2023
"""

##### Imports

from astropy.io import fits

import os
import saphires

# import numpy as np
# import pandas as pd

# import saphires as saph

# import glob
# import pickle
# import tqdm


##### Functions





##### Main script

def measure_radial_velocity( file_indices, header_df, config ):
    
    ### Loop through each of the files to measure RV for
    for i_file in file_indices:
        
        # Read in the spectra file
        file_name = os.path.join( config['paths']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) )
        file_in   = fits.open( file_name )
        
        # Continuum normalize the spectrum
        flux_cont_norm = file_in['extracted flux'].data / file_in['continuum'].data
        
        # Read in the target spectra!
        tar, tar_spec = saphires.io.read_vars( file_in['wavelength'].data, flux_cont_norm, header_df['file_token'].values[i_file], w_file = 'temp.ls', combine_all = False )

        # Prepare the spectra and run the broadening function!
        tar_spec = saphires.utils.prepare( tar, tar_spec, temp_spec, cr_trim = -0.2, oversample = 1 )

        try:
            tar_spec = saphires.bf.compute( tar, tar_spec, vel_width = config['radial_velocity']['bf_velocity_span'] )
        except RuntimeError:
            print( 'Issue with computing the broadening function!' )
            continue

        try:
            tar_spec = saphires.bf.analysis( tar, tar_spec, sb = 'sb1', R = config['radial_velocity']['bf_smooth_res'], single_plot = True, fit_trim = 20, text_out = False, prof = 'g' )

            os.rename( header_df['file_token'].values[i_file] + '.txt', night + '/bf/' + header_df['file_token'].values[i_file] + '/' + header_df['file_token'].values[i_file] + '.txt' )

            if designation[0] == '[':
                os.rename( '_allplots.pdf', night + '/bf/' + designation + '/' + designation[1:] + '_allplots.pdf' )
            else:
                os.rename( designation + '_allplots.pdf', night + '/bf/' + designation + '/' + designation + '_allplots.pdf' )

            saph.io.save( tar, tar_spec, night + '/bf/' + designation + '/' + designation + '_out' )

        except RuntimeError:
            print( 'Issue with fitting the broadening function!' )
            bad_objects.append( designation )
            os.rename( designation + '.txt', night + '/bf/' + designation + '/' + designation + '.txt' )

            saph.io.save( tar, tar_spec, night + '/bf/' + designation + '/' + designation + '_out' )
            continue


    return None