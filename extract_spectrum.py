""" Functions for extracting the echelle spectra.

Created by DMK on 10/22/2023
"""

##### Imports

import numpy as np

from astropy.io import fits
from scipy.optimize import curve_fit

import os

import tull_coude_utils

##### Functions

def get_order_image_block( full_image, order_width, order_trace ):
    
    # An empty array to hold the image block
    output_image_block = np.zeros( ( full_image.shape[1], order_width ) )
    
    # Go through each of the dispersion pixels and pull out the order image
    for pixel in range( full_image.shape[1] ):
        
        # Get the bottom and top edge of the order given the trace and cross dispersion order width
        bottom_edge = np.round( order_trace[pixel] ).astype( int ) - order_width // 2
        top_edge    = np.round( order_trace[pixel] ).astype( int ) + order_width // 2

        output_image_block[pixel] = full_image[bottom_edge:top_edge,pixel]
        
    return output_image_block

##### Main wrapper function

def extract_spectrum( file_indices, trace, header_df, config ):
    
    # If flagged in the config, reverse the order trace array
    if config['extraction']['reverse_traced_orders']:
        trace = trace[::-1]
        
    # The number of orders
    num_orders = trace.shape[0]
    
    # The number of dispersion pixels to extract
    num_pixels = trace.shape[1]

    ### Go through each of the files input to the function to extract
    for i_file in file_indices:
        
        ### Set up the frame to extract
        
        # The raw file name to read in and read it in
        file_name  = os.path.join( 'object_files', 'tullcoude_{}.fits'.format( header_df['file_token'].values[i_file] ) )
        input_file = fits.open( file_name )
        
        # Pull out the file's flux and error images
        full_image_flux = input_file[0].data        
        full_image_err  = input_file[1].data

        # Set up the output arrays for the extracted spectra
        extracted_flux = np.zeros( ( num_orders, num_pixels ) )
        extracted_err  = np.zeros( ( num_orders, num_pixels ) )
    
        ### Loop through the orders and do the actual extraction!
        for order in range( num_orders ):
            
            ### Pull out the flux and error image blocks around the trace
            order_image_flux = get_order_image_block( full_image_flux, config['extraction']['order_xdisp_width_extract'], trace[order] )
            order_image_err  = get_order_image_block( full_image_err, config['extraction']['order_xdisp_width_extract'], trace[order] )
        
        
    return None
