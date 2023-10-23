""" Library of miscellaneous functions used in the reduction pipeline for the Tull coude spectrograph
Last updated: DMK, 9/24/2023
"""

##### Imports

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.time import Time

import glob
import re

##### Functions

def make_header_manifest( header_manifest_file_name ):
    """
    
    Parameters
    ----------
    header_manifest_file_name : str
        File name for the output header manifest CSV file.

    Returns
    -------
    header_df : pandas DataFrame
        The data frame with the header information being written to a CSV file.

    """
    
    # Get the list of file names in the working directory
    file_names = np.sort( glob.glob( '*.fits' ) )
    
    # Set up the output data frame!
    header_df = pd.DataFrame( columns = [ 'file_name', 'file_token', 'image_type', 'object', 'ra', 'dec', 'exp_time', 'obs_jd', 'airmass', 'zenith_angle', 'gain', 'read_noise', 'em_flux'], 
                              index = range( file_names.size ) )

    # Go through each file and write header info to the output data frame
    for i_file, file_name in enumerate( file_names ):
        
        # Read in the header
        header = fits.getheader( file_name )
        
        # First write out the file_name
        header_df['file_name'][i_file] = file_name
        
        # Next write a file token
        file_token = re.sub( '[:-]', '', header['date-obs'] + 'T' + header['ut'].split('.')[0] )
        header_df['file_token'][i_file] = file_token
        
        # Image type
        header_df['image_type'][i_file] = header['imagetyp']
        
        # Object name
        header_df['object'][i_file] = header['object']
        
        # Coordinates
        if 'ra' in header:
            header_df['ra'][i_file] = header['ra']
        if 'dec' in header:
            header_df['dec'][i_file] = header['dec']
            
        # Exposure time
        header_df['exp_time'][i_file] = header['exptime']
        
        # The observation mid point
        header_df['obs_jd'][i_file] = Time( header['date-obs'] + 'T' + header['UT'], format = 'isot' ).jd
        
        # The airmass
        if 'airmass' in header:
            header_df['airmass'][i_file] = header['airmass']
            
        # The zenith angle
        header_df['zenith_angle'][i_file] = header['zd']
        
        # The gain and read noise -- check for if it is gain/rdnoise3 or gain/rdnoise2
        if 'gain3' in header:
            header_df['gain'][i_file], header_df['read_noise'][i_file] = header['gain3'], header['rdnoise3']
        else:
            header_df['gain'][i_file], header_df['read_noise'][i_file] = header['gain2'], header['rdnoise2']
            
        # E meter flux, if it is there!
        if 'emflux' in header:
            header_df['em_flux'][i_file] = header['emflux']

    # Write the header "manifest" out to a csv!
    header_df.to_csv( header_manifest_file_name, index = False )
    
    return header_df

def gaussian_1d( x_values, amplitude, mean, sigma, background ):
    """ Function to return 1D Gaussian.

    Parameters
    ----------
    x_values : array
        Array to evaluate the Gaussian at.
    amplitude : float
        Amplitude of the Gaussian (value at mean above the background).
    mean : float
        Mean of the Gaussian.
    sigma : float
        Standard deviation of the Gaussian.
    background : float
        Offset of the Gaussian.

    Returns
    -------
    y_values : array
        Array of Gaussian values at input x_values array.
    """
    
    y_values = amplitude * np.exp( - ( x_values - mean ) ** 2.0 / ( 2.0 * sigma ** 2.0 ) ) + background
    
    return y_values



