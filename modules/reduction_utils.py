""" Library of miscellaneous functions used in the reduction pipeline for the Tull coude spectrograph
Last updated: DMK, 9/24/2023
"""

##### Imports

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.time import Time
from scipy import stats

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

def polynomial_fit_sigma_reject( x_values, y_values, polynomial_degree, num_sigma_cut, num_iterations, y_limits = None, return_data = False ):
    """ Function to perform iterative polynomial fitting based on sigma rejection with the residuals
    
    Parameters
    ----------
    x_values : array
        Array of x values to fit
    y_values : array
        Array of y values to fit
    polynomial_degree : int
        Degree of the polynomial to fit
    num_sigma_cut : float
        The number of sigma beyond which residuals are rejected
    num_iterations : int
        The number of iterations for fitting
    y_limits : list, optional
        A list with the lower and upper limits for y data to fit. Default is None.
    return_data : bool, optional
        Flag for whether or not to return the x and y data arrays that are included in the final fit. Default is False.
    
    Returns
    -------
    poly_fit : array
        The best fit polynomial coefficients, the output of numpy polyfit
    x_values_to_fit : array, optional
        The x values included in the final fit that is output (after the iterative rejection). Returned if return_data is True.
    y_data_to_fit : array, optional
        The y values included in the final fit that is output (after the iterative rejection). Returned if return_data is True.
    """
    
    ### First -- get data to use in the first place: no nans and within parameter limits if they are given
    
    # Not nans
    not_nan = np.where( np.isfinite( y_values ) )[0]
    
    # Y data within limits if they are given!
    if y_limits is not None:
        within_limits = np.where( ( y_values > y_limits[0] ) & ( y_values < y_limits[1] ) )[0]
        
        # Intersect with the not nans to get data to use
        to_fit = np.intersect1d( not_nan, within_limits )

    # If no limits are given, only care about the not nans
    else:
        to_fit = not_nan
                
    x_values_to_fit = x_values[to_fit]
    y_values_to_fit = y_values[to_fit]
        
    ### Now do a first round fit!

    poly_fit = np.polyfit( x_values_to_fit, y_values_to_fit, polynomial_degree )
    
    ### Now loop through the number of iterations!
    
    for i_iter in range( num_iterations ):
        
        # Get the residuals
        fit_residuals = y_values_to_fit - np.polyval( poly_fit, x_values_to_fit )
        
        # The standard deviation (scaled MAD) of the residuals
        fit_residuals_stdev = stats.median_abs_deviation( fit_residuals, nan_policy = 'omit', scale = 'normal' )
                
        ### Another fit, where residuals are within the provide number of sigma for cutting
        
        to_fit = np.where( np.abs( fit_residuals ) < num_sigma_cut * fit_residuals_stdev )[0]
        
        # Re-define the arrays to fit
        x_values_to_fit = x_values_to_fit[to_fit]
        y_values_to_fit = y_values_to_fit[to_fit]
        
        poly_fit = np.polyfit( x_values_to_fit, y_values_to_fit, polynomial_degree )
    
    if return_data:
        return poly_fit, x_values_to_fit, y_values_to_fit
    else:
        return poly_fit

