""" Functions to continuum fit spectra with a spline.

Created by DMK on 11/13/2023
Last updated by DMK on 11/13/2023
"""

##### Imports
import numpy as np

from astropy.io import fits
from datetime import datetime
from scipy import interpolate, stats

import os

##### Functions

def continuum_fit_with_spline( x_values, y_values, x_knot_spacing, lower_sigma_reject, upper_sigma_reject, max_iter = 10 ):
    """ Function to fit a spline to a spectrum with sigma rejection. Different sigma rejection levels below and above the fit are allowed (i.e. get rid of absorption lines)
    It is assumed that no nans are present in the data provided.

    Parameters
    ----------
    x_values : array
        Array of x values to fit (independent variable, e.g. wavelength).
    y_values : array
        Array of y values to fit (dependent variable, e.g. flux).
    x_knot_spacing : float
        The knot spacing for the B spline.
    lower_sigma_reject : float
        The sigma level to reject points below the fit.
    upper_sigma_reject : float
        The sigma level to reject points above the fit.
    max_iter : int, optional
        The maximum number of sigma rejection iterations to perform. The default is 10.

    Returns
    -------
    spline_fit : tuple
        The tuple defining the best fit spline, as output by scipy.interpolate. Elements are (spline knot array, spline coefficient array, spline degree).
    """
    
    # Set the knots array. Keep in mind they must be interior knots, so start at the x value minimum + break space
    spline_knots = np.arange( x_values.min() + x_knot_spacing, x_values.max(), x_knot_spacing )
    
    # Set the values to fit array, will be modified in the rejection loop below
    x_values_to_fit = x_values.copy()
    y_values_to_fit = y_values.copy()

    # Loop for the maximum number of iterations, unless an iteration sooner results in no rejections
    for i_iter in range( max_iter ):

        # Get the b spline representation of the data
        spline_fit = interpolate.splrep( x_values_to_fit, y_values_to_fit, k = 3, t = spline_knots )

        # Calculate the residuals bewteen the data values and the spline fit
        residuals = y_values_to_fit - interpolate.splev( x_values_to_fit, spline_fit )

        # Get the standard deviation of the residuals, using the MAD
        residuals_mad = stats.median_abs_deviation( residuals, scale = 'normal' )

        # Keep points within the lower/upper sigma level provided
        within_mad = np.where( ( residuals < np.nanmedian( residuals ) + upper_sigma_reject * residuals_mad ) & ( residuals > np.nanmedian( residuals ) - lower_sigma_reject * residuals_mad ) )[0]

        # If no points are rejected, break out of the loop!
        if within_mad.size == y_values_to_fit.size:
            break

        # Re-define the x and y values to fit -- get rid of MAD rejected points
        x_values_to_fit = x_values_to_fit[within_mad]
        y_values_to_fit = y_values_to_fit[within_mad]
        
    return spline_fit

##### Main script function

def fit_spectra_continuum( file_indices, header_df, config ):
    """ Main script to run for fitting the continuum of science spectra. Fits the continuum with a spline and outputs it to a new FITS extension.

    Parameters
    ----------
    file_indices : array
        The file indices (in the header information file) of science observations with extracted spectra to measure RVs for.
    header_df : pandas DataFrame
        The compiled information from the file headers.
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    None.
    """
    
    for i_file in file_indices:
        
        # Read in the spectra file
        file_name = os.path.join( config['paths']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) )
        file_in   = fits.open( file_name )
        
        # Empty array for the continuum fit values
        continuum_values = np.full( ( file_in['extracted flux'].data.shape ), np.nan )
        
        # Go through each of the orders
        for order in range( continuum_values.shape[0] ):
            
            # Remove any flux nans
            not_nan = np.where( np.isfinite( file_in['extracted flux'].data[order] ) )[0]
            
            # Set the spline wavelength spacing based on the number of chunks to break the spectrum into (set in the config)
            wavelength_knot_spacing = round( np.ptp( file_in['wavelength'].data[order] ) / config['continuum_fit']['num_spectrum_chunks'] )

            # Do not let the knot spacing be smaller than the minimum value set in the config file
            if wavelength_knot_spacing < config['continuum_fit']['min_knot_spacing']:
                wavelength_knot_spacing = 15
            
            # Fit the continuum!
            continuum_spline = continuum_fit_with_spline( file_in['wavelength'].data[order,not_nan], file_in['extracted flux'].data[order,not_nan], wavelength_knot_spacing, config['continuum_fit']['lower_sigma_reject'], config['continuum_fit']['upper_sigma_reject'] )
            
            # Evaluate the spline
            continuum_values[order] = interpolate.splev( file_in['wavelength'].data[order], continuum_spline )
            
        ### Output the continuum fit
        
        # Build a new HDU list with the continuum extension added
        output_file = fits.HDUList( file_in[:4] + [ fits.ImageHDU( continuum_values, name = 'continuum' ) ] )
        
        # Add history to the primary header
        output_file[0].header['HISTORY'] = 'Continuum fit on {}'.format( datetime.strftime( datetime.now(), '%Y/%m/%d' ) )
        
        # Write out the file
        output_file.writeto( file_name, overwrite = True )
    
    return None
