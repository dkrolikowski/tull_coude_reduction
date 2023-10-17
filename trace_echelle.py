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

    # Empty array to hold the full trace
    full_trace = np.zeros( ( order_centers.size, flat_field_flux.shape[0] ), dtype = np.int )
    
    for pixel in range( 1, flat_field_flux.shape[0] + 1 ):
        
        if pixel == 1:
            prev_trace = order_centers
        else:
            prev_trace = full_trace[:,-pixel+1]
            
        flat_slice = flat_field_flux[:,-pixel+order_start_index+1]
        
        full_trace[:,-pixel] = recenter_order_trace( prev_trace, flat_slice, order_width_xdisp // 2 + 5 )
        
    plt.clf()
    plt.figure( figsize = ( 8, 6 ) )
    
    # Plot the image
    plt.imshow( np.log10( flat_field_flux ), cmap = 'gray', aspect = 'auto', interpolation = 'none' )
    
    for i_order in range( full_trace.shape[0] ):
        plt.plot( full_trace[i_order], '.', c = '#bf3465', ms = 1 )
    
    # Save the figure
    plt.savefig( os.path.join( config['paths']['reduction_dir'], 'plots', 'full_trace.pdf' ), bbox_inches = 'tight', pad_inches = 0.05 )
    plt.clf()
        
    return full_trace

def fit_full_trace( trace, flat_field_flux, trace_poly_degree, trace_poly_fit_start_index, config ):
    
    fit_trace = np.zeros( trace.shape )
    trace_poly_pars = np.zeros( ( trace.shape[0], trace_poly_degree + 1 ) )
    
    for i_order in range( trace.shape[0] ):
        trace_poly_pars[i_order] = np.polyfit( np.arange( trace_poly_fit_start_index, trace.shape[1] ), trace[i_order,trace_poly_fit_start_index:], trace_poly_degree )

    # Plot the initial trace fitting
    
    plt.clf()
    with PdfPages( os.path.join( config['paths']['reduction_dir'], 'plots', 'initial_fit_full_trace.pdf' ) ) as pdf:
        for y_range in [ ( 2048, 0 ), ( 2048, 1000 ), ( 950, 0 ) ]:
            plt.figure( figsize = ( 8, 6 ) )
            
            plt.imshow( np.log10( flat_field_flux ), aspect = 'auto', cmap = 'gray' )
            
            for i_order in range( trace.shape[0] ):
                plt.plot( trace[i_order], '.', c = '#dfa5e5', ms = 2 )
                
                plt.plot( np.arange( trace_poly_fit_start_index, trace[i_order].size ), trace[i_order,trace_poly_fit_start_index:], '.', c = '#50b29e', ms = 2 )
                
                plt.plot( np.polyval( trace_poly_pars[i_order], np.arange( trace.shape[1] ) ), '#bf3465', lw = 1 )
                
            plt.xlim( 0, 2048 )
            plt.ylim( y_range )
            
            pdf.savefig()
            plt.close()
            
    # Do a round of sigma clipping while hyper fitting the trace polynomial parameters!
    
    # hyper_poly_pars = np.zeros( ( 3, 4 ) )
    
    bad_orders = []
            
    # First go through the coefficients -- find orders that are "bad" across any of them
    for i_coeff in range( trace_poly_degree + 1 ):
        
        
        # Fit the trace polynomial coefficient across orders
        hyper_fit = np.polyfit( np.arange( trace_poly_pars.shape[0] ), trace_poly_pars[:,i_coeff], 2 )
        hyper_fit_vals = np.polyval( hyper_fit, np.arange( trace_poly_pars.shape[0] ) )
        
        # Sigma clip based on the absolute median deviation
        mad = np.nanmedian( np.abs( trace_poly_pars[:,i_coeff] - hyper_fit_vals ) )
        
        mask = ( np.abs( trace_poly_pars[:,i_coeff] - hyper_fit_vals ) <= 10 * mad )
        
        bad_orders.extend( np.where( ~mask )[0] )
        
    bad_orders  = np.unique( bad_orders )        
    good_orders = np.setdiff1d( np.arange( trace_poly_pars.shape[0] ), bad_orders )
        
    final_trace_poly_pars  = trace_poly_pars.copy()
    final_trace_hyper_pars = [ ]
    
    with PdfPages( os.path.join( config['paths']['reduction_dir'], 'plots', 'fit_trace_hyper_poly.pdf' ) ) as pdf:
        
        # Now fit only the "good" orders
        for i_coeff in range( trace_poly_degree + 1 ):
            
            fig, axs = plt.subplots( 2, 1, figsize = ( 10, 6 ), sharex = True, gridspec_kw = { 'height_ratios': [ 3, 1 ] } )

            # Second fit
            hyper_fit = np.polyfit( good_orders, trace_poly_pars[good_orders,i_coeff], 5 )
            hyper_fit_vals = np.polyval( hyper_fit, np.arange( trace_poly_pars.shape[0] ) )
            
            final_trace_hyper_pars.append( hyper_fit )
            
            axs[0].plot( np.arange( trace_poly_pars.shape[0] ), trace_poly_pars[:,i_coeff], 'o-' )
            axs[0].plot( good_orders, trace_poly_pars[good_orders,i_coeff], '*-' )

            axs[0].plot( np.arange( trace_poly_pars.shape[0] ), hyper_fit_vals, '-' )
            
            axs[1].plot( good_orders, np.polyval( hyper_fit, good_orders ) - trace_poly_pars[good_orders,i_coeff], 'o' )
            
            pdf.savefig()
            plt.close()
                       
            final_trace_poly_pars[bad_orders,i_coeff] = np.polyval( hyper_fit, bad_orders )
            
    plt.clf()
    with PdfPages( os.path.join( config['paths']['reduction_dir'], 'plots', 'final_fit_full_trace.pdf' ) ) as pdf:
        for y_range in [ ( 2048, 0 ), ( 2048, 1000 ), ( 950, 0 ) ]:
            plt.figure( figsize = ( 8, 6 ) )
            
            plt.imshow( np.log10( flat_field_flux ), aspect = 'auto', cmap = 'gray' )
            
            for i_order in range( trace.shape[0] ):
                plt.plot( trace[i_order], '.', c = '#dfa5e5', ms = 2 )
                
                plt.plot( np.arange( trace_poly_fit_start_index, trace[i_order].size ), trace[i_order,trace_poly_fit_start_index:], '.', c = '#874310', ms = 2 )
                
                plt.plot( np.polyval( trace_poly_pars[i_order], np.arange( trace.shape[1] ) ), '#bf3465', lw = 1 )
                if i_order in bad_orders:
                    plt.plot( np.polyval( final_trace_poly_pars[i_order], np.arange( trace.shape[1] ) ), '#50b29e', lw = 1, ls = '--' )

            plt.xlim( 0, 2048 )
            plt.ylim( y_range )
            
            pdf.savefig()
            plt.close()
            
    # Output the fit trace
    for i_order, fit_pars in enumerate( final_trace_poly_pars ):
        
        fit_trace[i_order] = np.polyval( fit_pars, np.arange( trace.shape[1] ) )
        
    # If the config says to extend the trace to 58 orders (and less than that are found) extend the hyper fit
    if fit_trace.shape[0] < config['trace']['number_of_orders']:
        
        fit_trace = np.zeros( ( config['trace']['number_of_orders'], trace.shape[1] ) )
        
        for i_order in range( config['trace']['number_of_orders'] ):
            
            if i_order < config['trace']['number_of_orders']:
                fit_trace[i_order] = np.polyval( final_trace_poly_pars[i_order], np.arange( trace.shape[1] ) )
            else:
                fit_pars = [ np.polyval( hyper_fit_pars, i_order ) for hyper_fit_pars in final_trace_hyper_pars ]
                fit_trace[i_order] = np.polyval( fit_pars, np.arange( trace.shape[1] ) )
        
    return fit_trace, bad_orders
    
##### Main wrapper script to find and fit the trace

def get_trace( flat_field, config ):
            
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
    fit_trace, orders_poorly_fit = fit_full_trace( full_trace, flat_field['flat flux'].data, config['trace']['trace_poly_degree'], config['trace']['trace_poly_fit_start_index'], config )
    
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
        


