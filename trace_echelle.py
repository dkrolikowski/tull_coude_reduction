""" Functions to find the echelle orders on the echellogram and fit their traces.

Created by DMK on 10/16/2023
Last updated by DMK on 10/16/2023
"""

##### Imports
import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
from scipy.stats import median_abs_deviation

##### Functions

def find_order_trace_center( flat_slice, method, order_width_xdisp = 20 ):
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
    
    for i_order, order in enumerate( order_trace_centers ):
        
        # Use a cut off of 0.7x the value of the flat slice at the order center to denote the "sides/edges"
        flat_cut_off = flat_slice[order] * 0.7
        
        # The padding to use on either side of the order center to try and catch the edges
        padding = order_width_xdisp // 2 + 5
        
        left_edges  = order - padding + np.where( flat_slice[order-padding:order] <= flat_cut_off )[-1]
        right_edges = order + np.where( flat_slice[order:order+padding+1] <= flat_cut_off )[-1]
        
        # If there are no indices in either of the left or right sides, just adopt the order center
        if len( left_edges ) == 0 or len( right_edges ) == 0:
            order_trace_centers[i_order] = order
        # Else adopt the halfway between the rightmost element of the left side and the leftmost element of the right side
        else:
            order_trace_centers[i_order] = ( left_edges[-1] + right_edges[0] ) // 2

    return order_trace_centers

##### Main wrapper script to find and fit the trace

def get_trace( config ):
    
    if config['trace']['do_step']:
        
        
        
    else:
        
    
    
    
    return trace





