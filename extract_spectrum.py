""" Functions for extracting the echelle spectra.

Created by DMK on 10/22/2023
"""

##### Imports

import numpy as np

from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.stats import median_abs_deviation

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

def optimal_extraction( flux_image, err_image, num_pixels, order_width, background_option ):
    
    ### Empty arrays to hold the output!
    flux_opt_extract = np.full( num_pixels, np.nan )
    err_opt_extract  = np.full( num_pixels, np.nan )
    
    ### First fit a background estimate. 
    # If the background value is fixed these are used to define the background. If the background value will be fit, this will provide the initial guess for the background

    # Fit a polynomial to the counts across the top of the order image. A 2nd order polynomial, iterating once for a 10 sigma rejection
    background_top_fit = tull_coude_utils.polynomial_fit_sigma_reject( np.arange( num_pixels ), np.mean( flux_image[:,:2], axis = 1 ), 2, 10, 1 )
    
    # Fit a polynomial to the counts across the bottom of the order image. A 2nd order polynomial, iterating once for a 10 sigma rejection
    background_bottom_fit = tull_coude_utils.polynomial_fit_sigma_reject( np.arange( num_pixels ), np.mean( flux_image[:,-2:], axis = 1 ), 2, 10, 1 )
    
    ### The first step: go through each dispersion pixel, and fit the cross dispersion flux slice with a Gaussian to get the spatial profile parameters over the order
        
    # Array to hold the fit parameters of the spatial Gaussian profile. If the background is fixed only 3 parameters, if the background is free then 4 parameters
    if background_option == 'fixed':
        spatial_profile_fit_pars = np.full( ( num_pixels, 3 ), np.nan )
    elif background_option == 'fit':
        spatial_profile_fit_pars = np.full( ( num_pixels, 4 ), np.nan )
    
    # Now loop through each of the dispersion pixels
    for pixel in range( num_pixels ):
    
        ### First pull out the order cross dispersion slice -- flux and error
        flux_slice = flux_image[pixel]
        err_slice  = err_image[pixel]
    
        ### Prepare to fit with a Gaussian
        
        # The x array is simply integer pixels from 0 to the width (fit x array will depend on nans)
        slice_xarr = np.arange( order_width )
    
        # Get rid of any nans! And if there are only 3 or fewer non-nans left -- don't fit just leave as nan!
        to_fit = np.where( np.isfinite( flux_slice ) )[0]
    
        if to_fit.size < 4:
            continue
    
        ### Now fit with a Gaussian
        
        # Evaluate the background estimate polynomial fits at this pixel. Either to use as a fixed background or inital guess for the background value
        background = 0.5 * ( np.polyval( background_top_fit, pixel ) + np.polyval( background_bottom_fit, pixel ) )
        
        # If the background is fixed to the values fit at the top of this function:
        if background_option == 'fixed':
    
            # Initial guess for the Gaussian parameters: Amp = avg of pix 7-9, Mean = pixel 8 (halfway), Sigma = 1.5 pixels
            p_guess = [ np.nanmean( flux_slice[7:10] ), 8, 1.5 ]
                
            # Create a bespoke Gaussian function with the fixed background
            gaussian_function_use = lambda x, a, m, s: tull_coude_utils.gaussian_1d( x, a, m, s, background )
            
        elif background_option == 'fit':
            
            # Initial guess for the Gaussian parameters: Amp = avg of pix 7-9, Mean = pixel 8 (halfway), Sigma = 1.5 pixels, Background = Value from fit
            p_guess = [ np.nanmean( flux_slice[7:10] ), 8, 1.5, background ]
                
            # Create a bespoke Gaussian function with the fixed background
            gaussian_function_use = tull_coude_utils.gaussian_1d

        # Fit with the Gaussian! Put inside a try/except in case the fit fails! And break out if it does
        try:
            gauss_fit_pars, _ = curve_fit( gaussian_function_use, slice_xarr[to_fit], flux_slice[to_fit], p0 = p_guess, sigma = err_slice[to_fit] )
        except:
            continue
            
        ### Evaluate the best fit Gaussian profile!
        
        # Get the bet fit values
        gauss_fit_vals = gaussian_function_use( slice_xarr, *gauss_fit_pars )
    
        # Calculate the residuals with the fit. In future might be used for rejecting bad pixels.
        gauss_fit_residuals     = flux_slice - gauss_fit_vals
        gauss_fit_residuals_mad = median_abs_deviation( gauss_fit_residuals, nan_policy = 'omit', scale = 'normal' )
    
        ### Make the normalized spatial Gaussian profile for optimal extraction
    
        # First just output the Gaussian fit parameters
        spatial_profile_fit_pars[pixel] = gauss_fit_pars.copy()
    
        # The amplitude is background subtracted and then normalized by the sum of the values. Background used depends on input option (either fixed value or from the fit)
        if background_option == 'fixed':
            normalized_profile_amp = ( gauss_fit_pars[0] - background ) / np.sum( gauss_fit_vals - background )
        elif background_option == 'fit':
            normalized_profile_amp = ( gauss_fit_pars[0] - gauss_fit_pars[-1] ) / np.sum( gauss_fit_vals - gauss_fit_pars[-1] )
    
        # Replace the output fit amplitude with the renormalized value
        spatial_profile_fit_pars[pixel,0] = normalized_profile_amp
        
    ### Now fit the spatial profile fit parameters cross dispersion pixels to enfore smooth behavior and corrected bad pixel slice fits
    # If the background is fixed, only the amplitude and sigma are fit
    # If the background is free in the Gaussian fit, it is also included here
    # The parameter indices are: Amplitude = 0, Sigma = 2, Background = 3. We fit with a 3rd degree polynomial.
    
    # The limits for the spatial profile parameter values to include in fits: amplitude between 0 and 1, sigma > 0, background has no limits
    spatial_profile_par_limits = { 0: [ 0, 1 ], 2: [ 0, np.inf ], 3: [ -np.inf, np.inf ] }

    # Empty arrays for the coefficients
    if background_option == 'fixed':
        fit_coeffs_spatial_profile_pars = { 0: np.zeros( 4 ), 2: np.zeros( 4 ) }
    elif background_option == 'fit':
        fit_coeffs_spatial_profile_pars = { 0: np.zeros( 4 ), 2: np.zeros( 4 ), 3: np.zeros( 4 ) }
            
    # Go through each of the parameters we are going to fit
    for i_par in fit_coeffs_spatial_profile_pars.keys():
    
        # Sigma reject polynomial fit! 3rd degree polynomial and 1 iteration of 5 sigma rejection.
        poly_fit = tull_coude_utils.polynomial_fit_sigma_reject( np.arange( num_pixels ), spatial_profile_fit_pars[:,i_par], 3, 5, 1, y_limits = spatial_profile_par_limits[i_par] )
    
        # Output the fit coefficients
        fit_coeffs_spatial_profile_pars[i_par] = poly_fit
        
    ### Now do the optimal extraction -- with the fits of spatial profile parameters vs dispersion pixel applied
    for pixel in range( num_pixels ):
    
        # Skip if the initial spatial Gaussian profile fit failed (is nan)
        if np.isnan( spatial_profile_fit_pars[pixel,0] ):
            continue
    
        ### First pull out the order cross dispersion slice information again
        flux_slice = flux_image[pixel]
        err_slice  = err_image[pixel]
        slice_xarr = np.arange( order_width )
    
        ### Now get the Gaussian profile
    
        # Initialize with the fit values
        slice_fit_pars = spatial_profile_fit_pars[pixel].copy()
    
        # Replace the parameters that need to be!
        for i_par, par_fit in fit_coeffs_spatial_profile_pars.items():
            slice_fit_pars[i_par] = np.polyval( par_fit, pixel )
    
        # Evaluate with the spatial profile fit parameters and re-normalize! The background provided here is always zero, since it was already subtracted for getting the amplitude
        gauss_fit_profile = tull_coude_utils.gaussian_1d( slice_xarr, *slice_fit_pars, 0 )
        gauss_fit_profile = gauss_fit_profile / np.nansum( gauss_fit_profile )
    
        ### Extract!
        
        # Get the background value -- either from the order image edge fits if fixed, or the spatial profile parameter fitting. In the latter case, set to 0 if negative
        if background_option == 'fixed':
            background = 0.5 * ( np.polyval( background_top_fit, pixel ) + np.polyval( background_bottom_fit, pixel ) )
        elif background_option == 'fit':
            background = slice_fit_pars[-1] if slice_fit_pars[-1] > 0 else 0

        # Subtract the background
        flux_slice_background_sub = flux_slice - background
    
        # Optimally extract
        flux_opt_extract[pixel] = np.nansum( flux_slice_background_sub * gauss_fit_profile / err_slice ** 2.0 ) / np.nansum( gauss_fit_profile ** 2 / err_slice ** 2.0 )
    
        # Revise slice error estimate
        err_slice_revise = np.sqrt( err_slice ** 2.0 + ( flux_opt_extract[pixel] * gauss_fit_profile + background ) )
    
        # The error on the extracted flux value
        err_opt_extract[pixel] = np.sqrt( np.nansum( gauss_fit_profile ) / np.nansum( gauss_fit_profile ** 2 / err_slice_revise ** 2.0 ) )

    return flux_opt_extract, err_opt_extract

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
        flux_full_image = input_file[0].data        
        err_full_image  = input_file[1].data

        # Set up the output arrays for the extracted spectra
        extracted_flux = np.zeros( ( num_orders, num_pixels ) )
        extracted_err  = np.zeros( ( num_orders, num_pixels ) )
    
        ### Loop through the orders and do the actual extraction!
        for order in range( num_orders ):
            
            ### Pull out the flux and error image blocks around the trace
            flux_order_image = get_order_image_block( flux_full_image, config['extraction']['order_xdisp_width_extract'], trace[order] )
            err_order_image  = get_order_image_block( err_full_image, config['extraction']['order_xdisp_width_extract'], trace[order] )
            
            if extraction_method == 'optimal':
                optimal_extraction( flux_order_image, err_order_image, config['extraction']['background_subtraction'] )
            elif extraction_method == 'sum':
                sum_extraction()
        
        
    return None
