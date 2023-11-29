""" Functions to find and trace the echelle orders, fit each order's trace with a polynomial, and write the trace to a file for use in the extraction module.

Created by DMK on 10/16/2023

Last updated by DMK on 10/18/2023
"""

##### Imports
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import find_peaks
from scipy.stats import median_abs_deviation

import os

##### Functions

def recenter_order_trace( centers, slice_values, padding ):
    """ Function to recenter a trace given an input spectral slice.
    The cross dispersion "edges" of the order are looked for within a set padding number of pixels from the input order trace center.
    As of 10/18/2023 the pipeline functions with this slice being the flat field, but in principle it could be any frame (such as a bright star)

    Parameters
    ----------
    centers : array
        Input array with initial order trace centers (shape: number of orders).
    slice_values : array
        The spectral cross-dispersion slice to use for re-centering.
    padding : int
        The padding radius (above/below) the order center in the cross dispersion direction to find the order edges.

    Returns
    -------
    re_centers : array
        The re-centered array of the order trace (shape: number of orders).
    """
        
    # Array to hold the output re-centered order traces
    re_centers = centers.copy()
    
    # Go through each of the orders that has been found
    for order, order_center in enumerate( centers ):
        
        # Use a cut off of 0.7x the value of the flat slice at the order center to denote the "sides/edges"
        flat_cut_off = slice_values[order_center] * 0.7
                
        # Get the indices of the slice on either side of the center within the padding that have flux values below the threshold cut off
        left_edges  = order_center - padding + np.where( slice_values[order_center-padding:order_center] <= flat_cut_off )[-1]
        right_edges = order_center + np.where( slice_values[order_center:order_center+padding+1] <= flat_cut_off )[-1]
        
        # If there are no such indices on either side of the order center, just adopt the input order center
        if len( left_edges ) == 0 or len( right_edges ) == 0:
            re_centers[order] = order_center
        # Otherwise, adopt the halfway value between the rightmost element of the left side and the leftmost element of the right side
        else:
            re_centers[order] = ( left_edges[-1] + right_edges[0] ) // 2

    return re_centers

def find_order_centers_along_slice( flat_slice, order_xdisp_width, plot_dir, method = 'peak_find' ):
    """ Function to find echelle order centers given a slice of the flat field along the cross dispersion axis. 
    The gradient of this flat slice is used to find the top edge of each order (where the flux steeply increases at the order edge).
    There are two methods available to use: one that uses a single threshold value for the flat slice gradient and one to use a scipy.signal peak finding algorithm.
    The former was from the original pipeline written with Aaron Rizzuto, although has been adjusted to use an estimate of the gradient's background noise level rather than a percentile cut.
    The latter is being added in this new version of the pipeline (as of October 2023).

    Parameters
    ----------
    flat_slice : array
        An array of flat field flux values along a cross dispersion axis slice (typically one pixel within the overscan region).
    order_xdisp_width : int
        A rough estimate of the cross dispersion width of an order. For our default Tull coude spectrograph setup this is roughly 20 pixels (order fully filled by flat field).
    plot_dir : str
        The path to where plots for the trace should be saved.
    method : str, optional
        One of either 'peak_find' or 'gradient_threshold'. The default is 'peak_find'

    Returns
    -------
    order_centers : array
        Array with the cross disperion pixel location of the order center (shape: number of orders found)
    """
    
    # Assert that the method is one of those allowed
    assert method in [ 'peak_find', 'gradient_threshold' ], "The order finding method must be either 'peak_find' or 'gradient_threshold'."

    # Calculate the gradient of the flat slice
    flat_slice_gradient = np.gradient( flat_slice )
    
    ### Check which method to use and apply the one that's selected
    if method == 'peak_find':
        
        # For this version we set values of the gradient below the median to the median: we don't care about the negative peaks (order's bottom edge) and it messes with the peak finding.
        flat_slice_gradient_use = flat_slice_gradient.copy()
        
        below_median_loc = np.where( flat_slice_gradient_use < np.nanmedian( flat_slice_gradient_use ) )[0]
        
        flat_slice_gradient_use[below_median_loc] = np.nanmedian( flat_slice_gradient_use )
        
        # Use scipy find_peaks with constraints on the distane between two order peaks (the input rough order width value) and the width of the gradient peak (2.01 to 4 pixels, from tests)
        order_centers, _ = find_peaks( flat_slice_gradient_use, distance = order_xdisp_width, width = [ 2.01, 4 ] )
        
        # Add half of the input order cross dispersion width so the estimates are for order center and not edge
        order_centers += order_xdisp_width // 2

    elif method == 'gradient_threshold':
        
        # Get the standard deviation at the right end of the flat field slice to estimate background noise. Use median absolute deviation for potential outliers
        gradient_stdev = median_abs_deviation( flat_slice_gradient[-250:-5], nan_policy = 'omit', scale = 'normal' )
        
        ## Now do the order finding!
        
        # Empty list to hold the order centers
        order_centers = []

        # "Empty" variable to hold the previous order center found
        last_order = 0
        
        # Loop through each pixel, start a few pixels in to avoid weird edge effects that can crop up
        for i_pix in range( 6, flat_slice.shape[0] ):
            
            # Check if the gradient of the flat slice is above the threshold: 20 * the estimated background noise
            if flat_slice_gradient[i_pix] > 20 * gradient_stdev:
                
                # Check that this pixel above the threshold is sufficiently far from the previously found order (input to function, but ~20 pixels), also make sure it isn't super far from last order (to prevent bad orders found at far end of slice)
                # If this is the first pixel above the threshold set to the first order
                if ( 100 > i_pix - last_order > order_xdisp_width ) or last_order == 0:
                    # Add the estimated half width of an order, since the left edge is being found
                    order_centers.append( i_pix + order_xdisp_width // 2 )
                    
        order_centers = np.array( order_centers )
                            
    ### Common to both methods -- recenter the found orders! More important for the gradient_threshold method but still worth doing just in case
    
    order_centers = recenter_order_trace( order_centers, flat_slice, order_xdisp_width // 2 + 5 )
                
    ### Plot the flat slice and the order centers that it found!
    
    plt.clf()
    plt.figure( figsize = ( 15, 6 ) )
    
    # Go through the orders that were found and plot their centers as vertical lines
    for order_center in order_centers:
        plt.axvline( x = order_center, ls = '-', lw = 1, c = '#dfa5e5' )

    # Plot the flat field slice used for finding the order centers
    plt.plot( flat_slice, '-', c = '#323232', lw = 0.75 )
    
    plt.gca().set_yscale( 'log' )
    
    # Labels and such
    
    plt.xlabel( 'Cross Dispersion Pixel' )
    plt.ylabel( 'Flat Field' )
    
    plt.title( '{} order centers found'.format( order_centers.size ) )
    
    plt.savefig( os.path.join( plot_dir, 'order_start_centers.pdf' ), bbox_inches = 'tight', pad_inches = 0.05 )
    plt.clf()

    return order_centers

def find_full_trace( flat_field_flux, order_centers, order_disp_start_index, order_xdisp_width, plot_dir ):
    """ Function to find the full trace from starting locations found at the edge of the echellogram.
    The orders are traced along the flat field, defined with roughly estimating the center of the flat field order 2D spectrum.

    Parameters
    ----------
    flat_field_flux : array
        The flat field flux image.
    order_centers : array
        Array with the cross disperion pixel location of the order centers from find_order_centers_along_slice (shape: number of orders found).
    order_disp_start_index : int
        The ~dispersion direction (x-axis) pixel index at which the starting order centers was found (typically to avoid the overscan region).
    order_xdisp_width : int
        A rough estimate of the cross dispersion width of an order in the flat field, defined in the config file.
    plot_dir : str
        The path to where plots for the trace should be saved.

    Returns
    -------
    full_trace : array
        The full trace found from the starting locations using the flat field (shape: number of orders, number of pixels).
    """
    
    # Empty array to hold the full trace -- shape: number of orders found to start, ~cross-dispersion (y-axis) pixels (no overscan)
    full_trace = np.zeros( ( order_centers.size, flat_field_flux.shape[0] ), dtype = np.int )
    
    # Set the previous trace location to be the input order centers
    prev_trace = order_centers
    
    # Loop through each pixel along the flat field image trace
    for pixel in range( 1, full_trace.shape[1] + 1 ):
        
        # Pull out the cross-dispersion direction slice of the flat field at this pixel
        flat_slice = flat_field_flux[:,-pixel+order_disp_start_index+1]
        
        # Recenter the trace using the previous pixel's center as a starting point, based on the flat field flux slice at the current pixel
        full_trace[:,-pixel] = recenter_order_trace( prev_trace, flat_slice, order_xdisp_width // 2 + 5 )
        
        # Set the previous trace location to be this pixel's trace
        prev_trace = full_trace[:,-pixel]
        
    ### Plot the full trace on top of the flat field image, with three different figure pages of different y-axis ranges
    
    plot_trace_on_image( flat_field_flux, full_trace, [ ( 2048, 0 ), ( 2048, 1000 ), ( 1050, 0 ) ], os.path.join( plot_dir, 'full_trace.pdf' ), title_str = 'Full Trace' )
    
    return full_trace

def fit_full_trace( trace, trace_poly_degree, trace_poly_fit_start_index, number_of_orders, flat_field_flux, plot_dir ):
    """ Function to take the full trace from find_full_trace and fit each order's trace with a polynomial.
    Then, the order trace polynomial coefficients themselves are fit as a function of order, and these hyper-fits are used to find "bad orders" with poorly defined traces.
    Poorly-fit orders have their trace polynomial fit coefficients replaced with coefficients from the hyper-fits.
    
    Parameters
    ----------
    trace : array
        The full trace found in find_full_trace (shape: number of orders, number of pixels)
    trace_poly_degree : int
        The polynomial degree to use for fitting each order's trace as defined in the config file.
    trace_poly_fit_start_index : int
        The starting pixel index for inclusion in fitting the trace as defined in the config file. This is not necessarily 0 because the trace could be poorly defined at the order ends.
    number_of_orders : int
        The config-defined number of orders to trace. If more orders are found, this function returns more orders than number_of_orders. If fewer orders are found, the trace is filled out using the polynomial fit coefficient hyper-fits.
    flat_field_flux : array
        The flat field flux image.
    plot_dir : str
        The path to where plots for the trace should be saved.

    Returns
    -------
    fit_trace : array
        The values of the trace resulting from the polynomial fits to each order's input full trace.
    bad_orders : array
        A list of orders with bad polynomial trace fits, which are replaced using hyper-fits to the polynomial coefficient of well-defined order traces.
    """
    
    ### Initial fit to the order traces
    
    # Set up empty arrays to hold the initial fit trace and the polynomial fit coefficients
    initial_fit_trace   = np.zeros( trace.shape )
    trace_poly_fit_pars = np.zeros( ( trace.shape[0], trace_poly_degree + 1 ) )
    
    # Go through each of the found orders and fit with a polynomial!
    for order in range( trace.shape[0] ):
        
        # Only fit the pixels from the config-defined start index onwards (at the left end of the echellogram the fit can get funky for bad orders)
        trace_poly_fit_pars[order] = np.polyfit( np.arange( trace_poly_fit_start_index, trace.shape[1] ), trace[order,trace_poly_fit_start_index:], trace_poly_degree )

        # Evaluate the fit
        initial_fit_trace[order] = np.polyval( trace_poly_fit_pars[order], np.arange( trace.shape[1] ) )
        
    # Plot the initial trace fit
    plot_file_name = os.path.join( plot_dir, 'initial_fit_trace.pdf' )
    plot_trace_on_image( flat_field_flux, trace, [ ( 2048, 0 ), ( 2048, 1000 ), ( 1050, 0 ) ], plot_file_name, title_str = 'Initial Trace Fit, using Pixels {}:'.format( trace_poly_fit_start_index ), trace_fit = initial_fit_trace )
    
    ### Fit the polynomial coefficients as a function of order, and sigma clip to find "bad orders" with poorly fit/defined traces to replace
    
    # Empty list to hold the bad orders
    bad_orders = []
            
    # Go through each of the polynomial coefficients and fit. Bad orders are those that are sigma clipped in ANY of the coefficients 
    for i_coeff in range( trace_poly_degree + 1 ):
        
        # Fit the trace polynomial coefficient across orders with a 2nd order polynomial (maybe make this a config defined variable at some point?)
        hyper_fit      = np.polyfit( np.arange( trace_poly_fit_pars.shape[0] ), trace_poly_fit_pars[:,i_coeff], 2 )
        hyper_fit_vals = np.polyval( hyper_fit, np.arange( trace_poly_fit_pars.shape[0] ) )
        
        ## Sigma clip based on the absolute median deviation
        mad  = np.nanmedian( np.abs( trace_poly_fit_pars[:,i_coeff] - hyper_fit_vals ) )
        
        # Orders with residuals greater than 10x the MAD are considered bad
        mask = ( np.abs( trace_poly_fit_pars[:,i_coeff] - hyper_fit_vals ) <= 10 * mad )
        
        bad_orders.extend( np.where( ~mask )[0] )
        
    # Get the unique list of bad orders (some may be bad in more than one coefficient) and define the good orders as all others
    bad_orders  = np.unique( bad_orders )        
    good_orders = np.setdiff1d( np.arange( trace_poly_fit_pars.shape[0] ), bad_orders )
    
    ### Now do a second round of polynomial coefficent fitting using only the good orders to get hyper-fits that will be used to define the trace for the bad orders
    
    # Arrays to hold the polynomial fit parameters to: the order traces and also the trace coefficient hyper-fits
    final_trace_poly_pars  = trace_poly_fit_pars.copy()
    final_trace_hyper_pars = np.zeros( ( trace_poly_degree + 1, 6 ) )
    
    # Wrap this in a multi-page plot with figures for each of the order trace polynomial coefficients and their hyper-fits
    with PdfPages( os.path.join( plot_dir, 'fit_trace_hyper_poly.pdf' ) ) as pdf:
        
        # Go through each of the coefficients and fit, using only the good orders
        for i_coeff in range( trace_poly_degree + 1 ):
            
            fig, axs = plt.subplots( 2, 1, figsize = ( 10, 6 ), sharex = True, gridspec_kw = { 'height_ratios': [ 3, 1 ] } )

            # Hyper fit on the good orders, using higher order polynomial than before (more accurate and now less prone to bad hyper-fit with the bad orders excluded)
            hyper_fit      = np.polyfit( good_orders, trace_poly_fit_pars[good_orders,i_coeff], 5 )
            hyper_fit_vals = np.polyval( hyper_fit, np.arange( trace_poly_fit_pars.shape[0] ) )
            
            # Plot the order trace polynomial coefficients, with the good orders highlighted as stars
            axs[0].plot( np.arange( trace_poly_fit_pars.shape[0] ), trace_poly_fit_pars[:,i_coeff], '.-', c = '#874310', ms = 7, lw = 0.85 )
            axs[0].plot( good_orders, trace_poly_fit_pars[good_orders,i_coeff], '*', c = '#dfa5e5', ms = 7, label = '"Good" Orders Included in Hyper Fit' )

            # Plot the hyper-fit to the coefficients
            axs[0].plot( np.arange( trace_poly_fit_pars.shape[0] ), hyper_fit_vals, '-', c = '#bf3465', lw = 1.5 )
            
            # Residuals on the bottom panel
            axs[1].plot( good_orders, np.polyval( hyper_fit, good_orders ) - trace_poly_fit_pars[good_orders,i_coeff], '.' )
            
            ## Labels and legend
            
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
    
    # Empty array to hold the final fit trace -- number of orders is the maximum between the config-defined number of orders and the number found using the starting method
    final_fit_trace = np.zeros( ( max( trace.shape[0], number_of_orders ), trace.shape[1] ) )
    
    # Go through each of the orders to fill out the fitted trace values
    for order in range( final_fit_trace.shape[0] ):
        
        # If the order is within the number of orders found in the full trace, adopt the polynomial fit to its trace
        if order < trace.shape[0]:
            fit_pars = final_trace_poly_pars[order]
        # If the order is beyond those found by the full trace, extend the hyper-fits to the trace polynomial coefficients
        else:
            fit_pars   = [ np.polyval( hyper_fit_pars, order ) for hyper_fit_pars in final_trace_hyper_pars ]
            # Append this to "bad" orders, for purpose of the below plot and output FITS header history
            bad_orders = bad_orders.append( order )
            
        final_fit_trace[order] = np.polyval( fit_pars, np.arange( trace.shape[1] ) )

    # for order, fit_pars in enumerate( final_trace_poly_pars ):
    #     final_fit_trace[order] = np.polyval( fit_pars, np.arange( trace.shape[1] ) )
    
    # Plot the final fit trace on the flat field
    plot_file_name = os.path.join( plot_dir, 'final_fit_trace.pdf' )
    plot_trace_on_image( flat_field_flux, trace, [ ( 2048, 0 ), ( 2048, 1000 ), ( 1050, 0 ) ], plot_file_name, title_str = 'Final Trace Fit, using Pixels {}:'.format( trace_poly_fit_start_index ), trace_fit = final_fit_trace, orders_to_highlight = bad_orders )
           
    return final_fit_trace, bad_orders

def plot_trace_on_image( image, trace, y_ranges, file_name, title_str = '', trace_fit = None, orders_to_highlight = [] ):
    """ Function to plot trace on top of an image. It can plot multiple figures in a multi-page PDF with different y ranges for zooming in on the echellogram.
    It can also take an array with the fit trace to overplot as a line for each order, and highlight certain orders (e.g. poorly fit orders) to be plotted as a different formatted line.

    Parameters
    ----------
    image : array
        The 2D image to plot behind the trace.
    trace : array
        The order trace points to overplot (shape: number of orders, number of pixels).
    y_ranges : list
        A list of tuples with the y-ranges to plot as individual figures in the multi-page PDF. The motivation is to allow for zoomed-in looks at the trace.
    file_name : str
        The file name and path for the output figure.
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
    
    # Set up the multi-page PDF (will only actually be multiple pages if the y_ranges list input has more than one entry).
    with PdfPages( file_name ) as pdf:
        
        # Loop through each of the y ranges input to the figure -- allowing for different pages at different zoom-in levels
        for y_range in y_ranges:
            plt.figure( figsize = ( 8, 6 ) )
            
            # Plot the background image, in log scale
            plt.imshow( np.log10( image ), aspect = 'auto', cmap = 'gray' )
            
            # Loop through each of the orders in the trace and over plot as points
            for order in range( trace.shape[0] ):
                plt.plot( trace[order], '.', c = '#dfa5e5', ms = 1 )
                
                # If a fitted trace to plot as lines is passed to the function
                if trace_fit is not None:
                    if order in orders_to_highlight: # Plot any orders to highlight with dashed lines and a different color
                        plt.plot( trace_fit[order], '--', c = '#50b29e', lw = 1.25, zorder = 12 )
                    else:
                        plt.plot( trace_fit[order], '-', c = '#bf3465', lw = 1, zorder = 12 )         
                
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
    
    # Pull out the path that points to the trace directory
    trace_dir = os.path.join( config['paths']['reduction_dir'], 'trace' )

    ### Find the starting location of the order traces
    
    # The slice of the flat field to pass to the order center finding function
    flat_field_slice = flat_field['flat flux'].data[:,config['trace']['order_start_index']]
    
    # Get the order centers
    order_start_centers = find_order_centers_along_slice( flat_field_slice, config['trace']['order_xdisp_trace_width'], trace_dir, config['trace']['order_center_method'] )
            
    # Get rid of the first order if it isn't fully on the detector with a bit of space on top (2/3 of the config defined cross dispersion order width)
    if order_start_centers[0] < config['trace']['order_xdisp_trace_width'] * 2 / 3:
        order_start_centers = order_start_centers[1:]
    
    ### Get the full trace across the full echellogram and fit each order with polynomials
    
    # Get the full trace using the flat field and the starting order centers
    full_trace = find_full_trace( flat_field['flat flux'].data, order_start_centers, config['trace']['order_start_index'], config['trace']['order_xdisp_trace_width'], trace_dir )
    
    # Fit the trace: the returned fit trace is from the evaluated polynomial fits to each order's trace. 
    # This also returns a list of orders with poor or missing order trace fits, so their trace is defined using fits to other orders' polynomial coefficients
    fit_trace, orders_poorly_fit = fit_full_trace( full_trace, config['trace']['trace_poly_degree'], config['trace']['trace_poly_fit_start_index'], config['trace']['number_of_orders'], flat_field['flat flux'].data, trace_dir )
    
    ### Output the trace
    
    # Only output up to the number of orders we want to extract (more orders may have been found!)
    full_trace = full_trace[:config['trace']['number_of_orders']]
    fit_trace  = fit_trace[:config['trace']['number_of_orders']]
    
    # Output the trace as a multi-extension FITS file, with both the full trace found along the flat field and the fitted trace values
    trace_hdu = fits.HDUList( [ fits.PrimaryHDU(), fits.ImageHDU( full_trace, name = 'data trace' ), fits.ImageHDU( fit_trace, name = 'fitted trace' ) ] )
        
    # Add various information to the header!
    trace_hdu[0].header['FILENAME'] = 'trace.fits'
    trace_hdu[0].header['NORDERS']  = ( config['trace']['number_of_orders'], 'Number of orders traced' )
    trace_hdu[0].header['POLYDEG']  = ( config['trace']['trace_poly_degree'], 'Polynomial degree fit to trace' )
    
    trace_hdu[0].header['HISTORY']  = 'Trace generated on {}.'.format( datetime.strftime( datetime.now(), '%Y/%m/%d' ) )
    trace_hdu[0].header['HISTORY']  = 'Values generated by polynomial fit to flat field trace.'
    trace_hdu[0].header['HISTORY']  = 'Initial order centers found using {} method'.format( config['trace']['order_center_method'] )
    trace_hdu[0].header['HISTORY']  = 'The full trace was fit with {} degree polynomial'.format( config['trace']['trace_poly_degree'] )
    trace_hdu[0].header['HISTORY']  = 'The trace was fit using x-axis (dispersion) pixels {} and up'.format(config['trace']['trace_poly_fit_start_index']  )
    trace_hdu[0].header['HISTORY']  = 'Orders with fit trace replaced using info from other order fit coeffs:'
    for order in np.intersect1d( np.arange( config['trace']['number_of_orders'] ), orders_poorly_fit ): 
        trace_hdu[0].header['HISTORY'] = 'Order {}'.format( order )
    
    trace_hdu.writeto( os.path.join( trace_dir, 'trace.fits' ), overwrite = True )
    
    return None

