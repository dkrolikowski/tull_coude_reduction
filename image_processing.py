"""
Library of functions to process images and write out the image files.

Created by DMK on 10/14/2023
Last updated by DMK on 10/14/2023
"""

##### Imports
import numpy as np

from astropy.io import fits
from astroscrappy import detect_cosmics
from datetime import datetime

import os
import pdb

##### Functions

##### Main wrapper function

def build_images( file_indices, super_bias, flat_field, bad_pixel_mask, header_df, config ):
    
    for i_file in file_indices:
        
        # The raw file name to read in
        file_name = header_df['file_name'].values[i_file]
        
        input_file = fits.open( file_name )
        
        # Assign the frame values, that can then be reassigned if cosmic subtraction is performed
        frame_values = input_file[0].data
                
        # Perform cosmic subtraction if flag is turned on in the config, and also only for object frames (not ThAr!)
        if config['image_process']['cosmic_subtract'] and input_file[0].header['imagetyp'] == 'object':
                        
            cosmics_mask, frame_values = detect_cosmics( input_file[0].data, sigclip = 5.0, sigfrac = 0.3, objlim = 5.0, niter = config['image_process']['cosmic_subtract_niter'], sepmed = True, satlevel = np.inf, 
                                                        cleantype = 'medmask', gain = header_df['gain'].values[i_file], readnoise = header_df['read_noise'].values[i_file] )
        
        # Calculate flux errors (in ADU)
        frame_errors = np.sqrt( frame_values + ( header_df['read_noise'].values[i_file] / header_df['gain'].values[i_file] ) ** 2.0 )
        
        # Now do basic image processing -- bias subtraction and flat fielding
        
        frame_values = frame_values - super_bias['bias flux'].data
        frame_values = frame_values / flat_field['flat flux'].data
        
        # Error
        frame_errors = np.sqrt( frame_errors ** 2.0 + super_bias['bias error'].data ** 2.0 + ( frame_values * flat_field['flat error'].data ) ** 2.0 ) / flat_field['flat flux'].data
        
        # Bad pixel mask
        bpm = np.where( bad_pixel_mask[0].data == 0 )
        
        # Set places where there is a bad pixel to the median of the frame flux, but with an extremely low S/N -- why not just nans?
        placeholder_bad_value = np.nanmedian( frame_values )
        placeholder_bad_error = np.nanmedian( frame_values ) / 1e-4
        
        frame_values[bpm] = placeholder_bad_value
        frame_errors[bpm] = placeholder_bad_error
        
        # And set any nans to the same
        nan_where = np.where( np.isnan( frame_values ) )
    
        frame_values[nan_where] = placeholder_bad_value
        frame_errors[nan_where] = placeholder_bad_error
        
        ### Set up the output image
        
        frame_value_hdu = fits.PrimaryHDU( frame_values, header = input_file[0].header )
        frame_error_hdu = fits.ImageHDU( frame_errors, name = 'error' )
        
        output_file = fits.HDUList( [ frame_value_hdu, frame_error_hdu ] )
        
        # Additions to the header
        
        output_file[0].header['FILENAME'] = 'tullcoude_{}.fits'.format( header_df['file_token'].values[i_file] )
        
        output_file[0].header['HISTORY'] = 'Image processed on {}'.format( datetime.strftime( datetime.now(), '%Y/%m/%d' ) )
        output_file[0].header['HISTORY'] = 'Bias subtracted and flat fielded'
        
        if config['image_process']['cosmic_subtract'] and input_file[0].header['imagetyp'] == 'object':
            output_file[0].header['HISTORY'] = 'Cosmic ray subtracted'
        
        output_file.writeto( os.path.join( config['paths']['reduction_dir'], 'object_files', 'tullcoude_{}.fits'.format( header_df['file_token'].values[i_file] ) ), overwrite = True )
        
    return None
