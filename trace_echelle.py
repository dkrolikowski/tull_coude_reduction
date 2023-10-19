""" Functions to find the echelle orders on the echellogram and fit their traces.

Created by DMK on 10/16/2023
Last updated by DMK on 10/16/2023
"""

##### Imports
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import find_peaks
from scipy.stats import median_abs_deviation

import os


##### Functions

def recenter_order_trace( centers, slice_values, padding ):
    """ Function to recenter a trace given an input spectral slice.
    The cross dispersion edges of the order are looked for within a set padding number of cross dispersion pixels from the input order trace center.

    Parameters
    ----------
    centers : array
        Input array with initial order trace centers (shape: number of orders).
    slice_values : array
        The spectral cross-dispersion slice to use for re-centering. Here we use a slice of the flat field to look for the top and bottom edges of the echelle orders.
    padding : int
        The padding to look away from the order center in the cross dispersion direction to find the order edges. This is slightly larger than the typical order width.

    Returns
    -------
    re_centers : array
        The re-centered array of the order trace (shape: number of orders).
    """
    
    ### Common to both methods -- recenter the found orders! More important for the gradient_threshold method but still worth doing just in case
    
    # Array to hold the output re-centered order traces
    re_centers = centers.copy()
    
    for i_order, order in enumerate( centers ):
        
        # Use a cut off of 0.7x the value of the flat slice at the order center to denote the "sides/edges"
        flat_cut_off = slice_values[order] * 0.7
                
        left_edges  = order - padding + np.where( slice_values[order-padding:order] <= flat_cut_off )[-1]
        right_edges = order + np.where( slice_values[order:order+padding+1] <= flat_cut_off )[-1]
                
        # If there are no indices in either of the left or right sides, just adopt the order center
        if len( left_edges ) == 0 or len( right_edges ) == 0:
            re_centers[i_order] = order
        # Else adopt the halfway between the rightmost element of the left side and the leftmost element of the right side
        else:
            re_centers[i_order] = ( left_edges[-1] + right_edges[0] ) // 2

    return re_centers

def find_order_centers_along_slice( flat_slice, method, order_width_xdisp, config ):
    """ Function to find the centers of orders given a slice of the flat field along the cross dispersion axis. The gradient of this flat slice is used to find the left edges of each order (where the flux steeply increases at the order edge).
    There are two methods currently coded for this: one that uses a single threshold value for the flat slice gradient, and one to use a scipy.signal peak finding algorithm.
    The former was from the original pipeline written with Aaron Rizzuto, although has been adjusted to use an estimate of the gradient's background noise level rather than a percentile cut.
    The latter is being added in this new version of the pipeline.

    Parameters
    ----------
    flat_slice : array
        An array of flat field flux values along a cross dispersion axis slice.
    method : str
        One of either 'peak_find' or 'gradient_threshold'.
    order_width_xdisp : int
        A rough estimate of the cross dispersion width of an order. For the Tull coude spectrograph this is roughly 20 pixels

    Returns
    -------
    order_trace_centers : array
        Array with the cross disperion pixel location of the order center.
    """
    
    # Assert that the method is one of those allowed
    assert method in [ 'gradient_threshold', 'peak_find' ], "The order trace center method must be either 'gradient threshold' or 'peak_find'."

    # Calculate the gradient of the flat slice
    flat_slice_gradient = np.gradient( flat_slice )
    
    ### Check which method to use and apply the one that's selected
    if method == 'gradient_threshold':
        
        # Get the standard deviation at the right end of the trace, to estimate background noise. Use median absolutey deviation in case it overlaps an order
        gradient_stdev = median_abs_deviation( flat_slice_gradient[-250:-5], nan_policy = 'omit', scale = 'normal' )
        
        ## Now do the order finding!
        order_trace_centers = []

        # "Empty" variable to hold the previous order center found
        last_order    = 0
        
        # Loop through each pixel, start a few pixels in to avoid some weird edge effects that can crop up
        for i_pix in range( 6, flat_slice.shape[0] ):
            
            # Check if the gradient of the flat slice is above the threshold: 20 * the estimated background noise
            if flat_slice_gradient[i_pix] > 20 * gradient_stdev:
                
                # Check that this pixel above the threshold is sufficiently far from the previously found order (input to function, but ~20 pixels), also make sure it isn't super far from last order (to prevent bad orders found at far end of slice)
                # If this is the first pixel above the threshold set to the first order
                if ( 100 > i_pix - last_order > order_width_xdisp ) or last_order == 0:
                    # Add the estimate half width of an order, since the left edge is being found
                    order_trace_centers.append( i_pix + order_width_xdisp // 2 )
                    
        order_trace_centers = np.array( order_trace_centers )
                    
    elif method == 'peak_find':
        
        # For this version, we set values of the gradient below the median to the median -- we don't care about the negative peaks (the right edge of the order) and it screws up the peak finding width parameter
        flat_slice_gradient_positive = flat_slice_gradient.copy()
        
        below_median_loc = np.where( flat_slice_gradient_positive < np.nanmedian( flat_slice_gradient_positive ) )[0]
        
        flat_slice_gradient_positive[below_median_loc] = np.nanmedian( flat_slice_gradient_positive )
        
        # Use scipy find peaks with constraints on the distane between two order peaks (the input rough order width value) and the width of the gradient peak (2.01 to 4 pixels, from tests)
        order_trace_centers, _ = find_peaks( flat_slice_gradient_positive, distance = order_width_xdisp, width = [ 2.01, 4 ] )
        
        # Add half of the input order cross dispersion width so the estimates are for order center and not edge
        order_trace_centers += order_width_xdisp // 2
        
    ### Common to both methods -- recenter the found orders! More important for the gradient_threshold method but still worth doing just in case
    
    order_trace_centers = recenter_order_trace( order_trace_centers, flat_slice, order_width_xdisp // 2 + 5 )
                
    ### Plot the flat slice and the order centers that it found!
    plt.clf()
    plt.figure( figsize = ( 15, 6 ) )
    
    for order_center in order_trace_centers:
        plt.axvline( x = order_center, ls = '-', lw = 1, c = '#dfa5e5' )

    plt.plot( flat_slice, '-', c = '#323232', lw = 0.75 )
    
    plt.gca().set_yscale( 'log' )
    
    plt.xlabel( 'Cross Dispersion Pixel' )
    plt.ylabel( 'Flat Field' )
    
    plt.title( '{} order centers found'.format( order_trace_centers.size ) )
    
    plt.savefig( os.path.join( config['paths']['reduction_dir'], 'plots', 'order_start_centers.pdf' ), bbox_inches = 'tight', pad_inches = 0.05 )
    plt.clf()

    return order_trace_centers

def find_full_trace( flat_field_flux, order_centers, order_start_index, order_width_xdisp, config ):
    """ Function to find the full trace from starting locations found at the edge of the echellogram.
    The orders are traced along the flat field, defined with roughly estimating the center of the flat field order 2D spectrum.

    Parameters
    ----------
    flat_field_flux : array
        The full 2D flat field flux array.
    order_centers : array
        Array with the cross disperion pixel location of the order centers from find_order_centers_along_slice (shape: number of orders found).
    order_start_index : int
        The index at which the starting order centers was found (typically to avoid the overscan region).
    order_width_xdisp : int
        A rough estimate of the cross dispersion width of an order in the flat field, defined in the config file.
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    full_trace : array
        The full trace found from the starting locations using the flat field (shape: number of orders found, number of pixels).
    """
    
    ### Get the full trace along the flat field

    # Empty array to hold the full trace -- shape: number of orders found to start, "y"-axis pixels (no overscan)
    full_trace = np.zeros( ( order_centers.size, flat_field_flux.shape[0] ), dtype = np.int )
    
    # Set the previous trace location to be to order centers
    prev_trace = order_centers
    
    # Loop through each pixel along the trace
    for pixel in range( 1, full_trace.shape[1] + 1 ):
        
        # Pull out the cross-dispersion direction slice of the flat field at this pixel
        flat_slice = flat_field_flux[:,-pixel+order_start_index+1]
        
        # Recenter the trace using the previous pixel's center as a starting point, based on the flat field flux slice at the current pixel
        full_trace[:,-pixel] = recenter_order_trace( prev_trace, flat_slice, order_width_xdisp // 2 + 5 )
        
        # Set the previous trace location to be this pixel's trace
        prev_trace = full_trace[:,-pixel]
        
    ### Plot the full trace on top of the flat field image
    
    plot_trace_on_image( flat_field_flux, full_trace, [ ( 2048, 0 ), ( 2048, 1000 ), ( 1050, 0 ) ], os.path.join( config['paths']['reduction_dir'], 'plots', 'full_trace.pdf' ), title_str = 'Full Trace' )
    
    return full_trace

def fit_full_trace( trace, trace_poly_degree, trace_poly_fit_start_index, flat_field_flux, config ):
    """ Function to take the full trace from find_full_trace and fit each order with a polynomial.
    Then, the polynomial coefficients are fit as a function of order, and these hyper-fits are chosen to find "bad orders" with poorly defined traces.
    Bad orders have their trace polynomial fits replaced with coefficients from these coefficient hyper-fits.
    
    Parameters
    ----------
    trace : array
        The full trace found in find_full_trace.
    trace_poly_degree : int
        The polynomial degree to use for fitting each order's trace as defined in the config file.
    trace_poly_fit_start_index : int
        The starting pixel index for inclusion in fitting the trace as defined in the config file. This would not be 0 if the end of the trace is poorly defined from find_full_trace.
    flat_field_flux : array
        The full 2D flat field flux array.
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    fit_trace : array
        The values of the trace resulting from the polynomial fits to each order's input full trace.
    bad_orders : array
        A list of orders with bad polynomial trace fits, which are replaced using hyper-fits to the polynomial coefficient of well-defined order traces.
    """
    
    ### Initial fit to the order traces
    
    # Set up empty arrays to hold the initial fit trace and the polynomial fit coefficients
    initial_fit_trace = np.zeros( trace.shape )
    trace_poly_pars   = np.zeros( ( trace.shape[0], trace_poly_degree + 1 ) )
    
    # Go through each of the found orders and fit with a polynomial!
    for i_order in range( trace.shape[0] ):
        
        # Only fit the pixels from the config-defined start index onwards (at the left end of the echellogram the fit can get funky for bad orders)
        trace_poly_pars[i_order] = np.polyfit( np.arange( trace_poly_fit_start_index, trace.shape[1] ), trace[i_order,trace_poly_fit_start_index:], trace_poly_degree )

        # Evaluate the fit
        initial_fit_trace[i_order] = np.polyval( trace_poly_pars[i_order], np.arange( trace.shape[1] ) )
        
    # Plot the initial trace fit
    plot_file_name = os.path.join( config['paths']['reduction_dir'], 'plots', 'initial_fit_trace.pdf' )
    plot_trace_on_image( flat_field_flux, trace, [ ( 2048, 0 ), ( 2048, 1000 ), ( 1050, 0 ) ], plot_file_name, title_str = 'Initial Trace Fit, using Pixels {}:'.format( trace_poly_fit_start_index ), trace_fit = initial_fit_trace )
    
    ### Fit the polynomial coefficients as a function of order, and sigma clip to find "bad orders" with poorly fit/defined traces to replace
    
    # Empty list to hold the bad orders
    bad_orders = []
            
    # Go through each of the polynomial coefficients and fit. Bad orders are those that are sigma clipped in ANY of the coefficients 
    for i_coeff in range( trace_poly_degree + 1 ):
        
        # Fit the trace polynomial coefficient across orders with a 2nd order polynomial (maybe make this a config defined variable at some point?)
        hyper_fit      = np.polyfit( np.arange( trace_poly_pars.shape[0] ), trace_poly_pars[:,i_coeff], 2 )
        hyper_fit_vals = np.polyval( hyper_fit, np.arange( trace_poly_pars.shape[0] ) )
        
        ## Sigma clip based on the absolute median deviation
        mad  = np.nanmedian( np.abs( trace_poly_pars[:,i_coeff] - hyper_fit_vals ) )
        
        # Orders with residuals greater than 10x the MAD are considered bad
        mask = ( np.abs( trace_poly_pars[:,i_coeff] - hyper_fit_vals ) <= 10 * mad )
        
        bad_orders.extend( np.where( ~mask )[0] )
        
    # Get the unique list of bad orders (some may be bad in more than one coefficient) and define the good orders as all others
    bad_orders  = np.unique( bad_orders )        
    good_orders = np.setdiff1d( np.arange( trace_poly_pars.shape[0] ), bad_orders )
    
    ### Now do a second round of polynomial coefficent fitting using only the good orders, to get hyper-fits that will be used to define the trace for the bad orders
    
    final_trace_poly_pars  = trace_poly_pars.copy()
    final_trace_hyper_pars = np.zeros( ( trace_poly_degree + 1, 6 ) )
    
    # Wrap this in a multi-page plot with figures for each of the order trace polynomial coefficients and their hyper-fits
    with PdfPages( os.path.join( config['paths']['reduction_dir'], 'plots', 'fit_trace_hyper_poly.pdf' ) ) as pdf:
        
        # Go through each of the coefficients and fit, using only the good orders
        for i_coeff in range( trace_poly_degree + 1 ):
            
            fig, axs = plt.subplots( 2, 1, figsize = ( 10, 6 ), sharex = True, gridspec_kw = { 'height_ratios': [ 3, 1 ] } )

            # Hyper fit on the good orders, using higher order polynomial than before (more accurate and now less prone to bad hyper-fit with the bad orders excluded)
            hyper_fit      = np.polyfit( good_orders, trace_poly_pars[good_orders,i_coeff], 5 )
            hyper_fit_vals = np.polyval( hyper_fit, np.arange( trace_poly_pars.shape[0] ) )
            
            
            # Plot the order trace polynomial coefficients, with the good orders highlight as stars
            axs[0].plot( np.arange( trace_poly_pars.shape[0] ), trace_poly_pars[:,i_coeff], '.-', c = '#874310', ms = 7, lw = 0.85 )
            axs[0].plot( good_orders, trace_poly_pars[good_orders,i_coeff], '*', c = '#dfa5e5', ms = 7, label = '"Good" Orders Included in Hyper Fit' )

            # Plot the hyper-fit to the coefficients
            axs[0].plot( np.arange( trace_poly_pars.shape[0] ), hyper_fit_vals, '-', c = '#bf3465', lw = 1.5 )
            
            axs[1].plot( good_orders, np.polyval( hyper_fit, good_orders ) - trace_poly_pars[good_orders,i_coeff], '.' )
            
            # Labels and legend
            
            axs[0].set_ylabel( 'Trace Polynomial Coefficient $c_{}$'.format( trace_poly_degree - i_coeff ) )
            axs[0].legend( fontsize = 'small' )
            
            axs[1].set_xlabel( 'Echelle Order Index' )
            axs[1].set_ylabel( 'Residuals' )
            
            pdf.savefig()
            plt.close()
            
            # Output the coefficient hyper-fit parameters
            final_trace_hyper_pars[i_coeff] = hyper_fit
            
            # Change the order trace polynomial coefficients to those from the hyper-fit for any bad orders
            if bad_orders.size > 0:
                final_trace_poly_pars[bad_orders,i_coeff] = np.polyval( hyper_fit, bad_orders )
            
    ### Get the final fit order trace values to output
    
    # Empty array to hold the final fit trace
    final_fit_trace = np.zeros( trace.shape )
    
    for i_order, fit_pars in enumerate( final_trace_poly_pars ):
        final_fit_trace[i_order] = np.polyval( fit_pars, np.arange( trace.shape[1] ) )
            
    # Plot the final fit trace on the flat field
    plot_file_name = os.path.join( config['paths']['reduction_dir'], 'plots', 'final_fit_trace.pdf' )
    plot_trace_on_image( flat_field_flux, trace, [ ( 2048, 0 ), ( 2048, 1000 ), ( 1050, 0 ) ], plot_file_name, title_str = 'Final Trace Fit, using Pixels {}:'.format( trace_poly_fit_start_index ), trace_fit = final_fit_trace, orders_to_highlight = bad_orders )

    # # If the config says to extend the trace to 58 orders (and less than that are found) extend the hyper fit
    # if fit_trace.shape[0] < config['trace']['number_of_orders']:
        
    #     fit_trace = np.zeros( ( config['trace']['number_of_orders'], trace.shape[1] ) )
        
    #     for i_order in range( config['trace']['number_of_orders'] ):
            
    #         if i_order < config['trace']['number_of_orders']:
    #             fit_trace[i_order] = np.polyval( final_trace_poly_pars[i_order], np.arange( trace.shape[1] ) )
    #         else:
    #             fit_pars = [ np.polyval( hyper_fit_pars, i_order ) for hyper_fit_pars in final_trace_hyper_pars ]
    #             fit_trace[i_order] = np.polyval( fit_pars, np.arange( trace.shape[1] ) )
        
    return final_fit_trace, bad_orders

def plot_trace_on_image( image, trace, y_ranges, file_name, title_str = '', trace_fit = None, orders_to_highlight = [] ):
    """ Function to plot a trace on top of an image. It can plot plot multiple images in a single multi-page PDF with different y ranges for zooming in on different echellogram areas.
    It can also take an array with the fit trace to overplot as a line for each order, and highlight certain orders (e.g. poorly fit orders) to be plotted as a different formatted line.

    Parameters
    ----------
    image : array
        The 2D image to plot in the background behind the trace.
    trace : array
        The order trace points to overplot (shape: number of orders, number of pixels).
    y_ranges : list
        A list of tuples with the y-ranges to plot as individual figures in the multi-page PDF. The motivation is to allow for zoomed-in looks at the trace.
    file_name : str
        The file name for the output figure.
    title_str : str, optional
        An additional string for the plot title to append to the title with information about the number of orders traced. The default is ''.
    trace_fit : array, optional
        An array of order trace fit values to overplot as a line for each order if provided. The default is None.
    orders_to_highlight : list, optional
        A list of orders to highlight in the fit trace (if provided) with dashed lines rather than solid lines. The default is [].

    Returns
    -------
    None.
    """
    
    # Set up the multi-page PDF (only will be multiple pages if the y_ranges list input has more than one entry).
    with PdfPages( file_name ) as pdf:
        
        # Loop through each of the y ranges input to the figure -- allowing for different pages at different zoom-in levels
        for y_range in y_ranges:
            plt.figure( figsize = ( 8, 6 ) )
            
            # Plot the background image, in log scale
            plt.imshow( np.log10( image ), aspect = 'auto', cmap = 'gray' )
            
            # Loop through each of the orders in the trace and over plot as points
            for i_order in range( trace.shape[0] ):
                plt.plot( trace[i_order], '.', c = '#dfa5e5', ms = 1 )
                
                if trace_fit is not None:
                    if i_order in orders_to_highlight:
                        plt.plot( trace_fit[i_order], '--', c = '#50b29e', lw = 1.25, zorder = 12 )
                    else:
                        plt.plot( trace_fit[i_order], '-', c = '#bf3465', lw = 1, zorder = 12 )         
                
            # Set the limits -- x limit is full image, y limit is input y range                
            plt.xlim( 0, image.shape[0] )
            plt.ylim( y_range )
            
            plt.title( 'Number of orders traced: {}, {}'.format( trace.shape[0], title_str ) )
            
            pdf.savefig( bbox_inches = 'tight', pad_inches = 0.05 )
            plt.close()

    return None
    
##### Main wrapper script to find and fit the trace

def get_trace( flat_field, config ):
    """ Main wrapper function to find and fit the full echellogram order traces, and output the fitted trace file for extracting spectra.
    It does not return the trace location array, but writes it to a FITS file.

    Parameters
    ----------
    flat_field : HDUList
        The FITS file containing the flat field.
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    None.
    """
            
    ### Find the starting location of the order traces
    
    # The slice of the flat field to pass to the order center finding function
    flat_field_slice = flat_field['flat flux'].data[:,config['trace']['order_start_index']]
    
    # Get the order centers
    order_start_centers = find_order_centers_along_slice( flat_field_slice, config['trace']['order_center_method'], config['instrument']['order_xdisp_width'], config )
            
    # Get rid of the first order if it isn't fully on the detector
    if order_start_centers[0] < 15:
        order_start_centers = order_start_centers[1:]
    
    # Get the full trace and fit it with polynomials per order
    full_trace = find_full_trace( flat_field['flat flux'].data, order_start_centers, config['trace']['order_start_index'], config['instrument']['order_xdisp_width'], config )
        
    # Fit the trace
    fit_trace, orders_poorly_fit = fit_full_trace( full_trace, config['trace']['trace_poly_degree'], config['trace']['trace_poly_fit_start_index'], flat_field['flat flux'].data, config )
    
    # Only output up to the number of orders we want to extract
    fit_trace = fit_trace[:config['trace']['number_of_orders']]
    
    # Ouptut the trace file
    trace_hdu = fits.PrimaryHDU( fit_trace )
    
    # Add to the header
    trace_hdu.header['FILENAME'] = 'fitted_trace.fits'
    trace_hdu.header['NORDERS']  = ( config['trace']['number_of_orders'], 'Number of orders traced' )
    trace_hdu.header['POLYDEG']  = ( config['trace']['trace_poly_degree'], 'Polynomial degree fit to trace' )
    
    trace_hdu.header['HISTORY']  = 'Trace generated on {}.'.format( datetime.strftime( datetime.now(), '%Y/%m/%d' ) )
    trace_hdu.header['HISTORY']  = 'Values generated by polynomial fit to flat field trace.'
    trace_hdu.header['HISTORY']  = 'Initial order centers found using {} method'.format( config['trace']['order_center_method'] )
    trace_hdu.header['HISTORY']  = 'The full trace was fit with {} degree polynomial'.format( config['trace']['trace_poly_degree'] )
    trace_hdu.header['HISTORY']  = 'The trace was fit using x-axis (dispersion) pixels {} and up'.format(config['trace']['trace_poly_fit_start_index']  )
    trace_hdu.header['HISTORY']  = 'Poorly fit orders replaced using hyper-fits to the poly coeffs are:'
    for order in np.intersect1d( np.arange( config['trace']['number_of_orders'] ), orders_poorly_fit ): 
        trace_hdu.header['HISTORY'] = 'Order {}'.format( order )
    
    
    trace_hdu.writeto( os.path.join( config['paths']['reduction_dir'], 'trace', 'trace.fits' ), overwrite = True )
    
    return None
        


