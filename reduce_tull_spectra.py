#!/usr/bin/env python

import numpy as np

from astropy.io import fits

import ccd_calibrations
import image_processing
import trace_echelle
import extract_spectrum
import wavelength_solve_and_calibrate
import continuum_fit

import tull_coude_utils
import yaml
import os

##### Read in the config file defining this reduction run

config_file = yaml.safe_load( open( 'test_config.yml', 'r' ) )

print( 'PIPELINE START: Running reduction pipeline for night at data path: {}'.format( config_file['paths']['working_dir'] ) )

##### Set up directories for pipeline output

# cd into the working directory
os.chdir( config_file['paths']['working_dir'] )

# Make the directory tree (will move this to be defined by the config file in the future)
directories = [ 'cals', 'trace', 'object_files', 'spectrum_files', 'wavecal' ]
for dir_name in directories:
    os.makedirs( os.path.join( config_file['paths']['reduction_dir'], dir_name ), exist_ok = True )

##### Pull header information from files and output

header_info = tull_coude_utils.make_header_manifest( config_file['paths']['header_manifest_file_name'] )

lower_case_object_names = np.array( [ object_name.lower() for object_name in header_info['object'].values ] )

# Build CCD calibrations!

# Get the file indices for the bias and flat frames. Done in the main script rather than function in case this requires some changes to the file sorting.
bias_frame_indices = np.where( header_info['image_type'].values == 'zero' )[0]
flat_frame_indices = np.where( ( header_info['image_type'].values == 'flat' ) & ( header_info['object'].values != 'FF integration time test' ) )[0]

if config_file['calibrations']['do_step']:
    print( 'MODULE: Making CCD calibration files.' )
    
    ccd_calibrations.build_calibrations( header_info, bias_frame_indices, flat_frame_indices, config_file )

# Generate the image files, after bias and flat processing

# Get the file indices for the object frames to process -- ThAr reference and science frames

thar_names = [ 'thar', 'a' ] # Lower case!

thar_frame_indices = np.where( np.logical_and( ( header_info['image_type'].values == 'comp' ), np.any( [ lower_case_object_names == name for name in thar_names ], axis = 0 ) ) )[0]

not_object_names = [ 'test' ]

object_frame_indices = np.where( np.logical_and( header_info['image_type'].values == 'object', np.all( [ lower_case_object_names != name for name in not_object_names ], axis = 0 ) ) )[0]

frames_to_extract = np.sort( np.concatenate( [ thar_frame_indices, object_frame_indices ] ) )

# If the config file indicates to generate the processed object images
if config_file['image_process']['do_step']:
    print( 'MODULE: Image processing.' )
    
    # Read in the CCD calibration files necessary for image processing if we are going to image process!
    super_bias     = fits.open( os.path.join( config_file['paths']['reduction_dir'], 'cals', 'super_bias.fits' ) )
    flat_field     = fits.open( os.path.join( config_file['paths']['reduction_dir'], 'cals', 'flat_field.fits' ) )
    bad_pixel_mask = fits.open( os.path.join( config_file['paths']['reduction_dir'], 'cals', 'bad_pixel_mask.fits' ) )
    
    image_processing.build_images( frames_to_extract, super_bias, flat_field, bad_pixel_mask, header_info, config_file )

#### Trace

# If the config file indicates to make the trace
if config_file['trace']['do_step']:
    print( 'MODULE: Tracing echelle orders.' )
    
    # Read in the flat field for getting the trace if we need to
    flat_field = fits.open( os.path.join( config_file['paths']['reduction_dir'], 'cals', 'flat_field.fits' ) )
    
    trace_echelle.get_trace( flat_field, config_file )

#### Extract!

if config_file['extraction']['do_step']:
    print( 'MODULE: Extracting 1D spectra.' )
    
    trace_file = fits.open( os.path.join( config_file['paths']['reduction_dir'], 'trace', 'trace.fits' ) )
        
    not_object_names = [ 'test', 'solar', 'sol port', 'solar port', 'solport', 'solarport', 'solar port halpha' ]

    object_frame_indices_to_extract = np.where( np.logical_and( header_info['image_type'].values == 'object', np.all( [ lower_case_object_names != name for name in not_object_names ], axis = 0 ) ) )[0]
    
    extract_spectrum.extract_spectrum( object_frame_indices_to_extract, trace_file[2].data, header_info, config_file['extraction']['science_extract_type'], config_file['extraction']['science_background_subtract'], config_file )

    thar_names = [ 'thar', 'a' ] # Lower case!
    
    thar_frame_indices_to_extract = np.where( np.logical_and( ( header_info['image_type'].values == 'comp' ), np.any( [ lower_case_object_names == name for name in thar_names ], axis = 0 ) ) )[0]

    extract_spectrum.extract_spectrum( thar_frame_indices_to_extract, trace_file[2].data, header_info, config_file['extraction']['lamp_extract_type'], config_file['extraction']['lamp_background_subtract'], config_file )
    
### Wave cal!

if config_file['wavecal']['do_step']:
    print( 'MODULE: Wavelength calibration.' )
    
    thar_names = [ 'thar', 'a' ] # Lower case!
    
    thar_frame_indices_to_wavesol = np.where( np.logical_and( ( header_info['image_type'].values == 'comp' ), np.any( [ lower_case_object_names == name for name in thar_names ], axis = 0 ) ) )[0]
    thar_frame_indices_to_wavesol = np.intersect1d( thar_frame_indices_to_wavesol, np.where( header_info['exp_time'].values >= config_file['wavecal']['min_arc_exp_time'] )[0] )
    
    wavelength_solve_and_calibrate.wavelength_solution( thar_frame_indices_to_wavesol, header_info, config_file )

    frame_indices_to_wavecal = [ i for i in range( header_info.shape[0] ) if os.path.exists( os.path.join( config_file['paths']['reduction_dir'], 'spectrum_files/tullcoude_{}_spectrum.fits'.format( header_info['file_token'].values[i] ) ) ) ]
    frame_indices_to_wavecal = np.setdiff1d( frame_indices_to_wavecal, thar_frame_indices_to_wavesol )
        
    wavelength_solve_and_calibrate.wavelength_calibrate( frame_indices_to_wavecal, thar_frame_indices_to_wavesol, header_info, config_file )
    
### Continuum fit

if config_file['continuum_fit']['do_step']:
    print( 'MODULE: Fitting continuum' )
    
    not_object_names = [ 'test', 'solar', 'sol port', 'solar port', 'solport', 'solarport', 'solar port halpha' ]
    object_frame_indices_to_extract = np.where( np.logical_and( header_info['image_type'].values == 'object', np.all( [ lower_case_object_names != name for name in not_object_names ], axis = 0 ) ) )[0]
    
    continuum_fit.fit_spectra_continuum( object_frame_indices_to_extract, header_info, config_file )
        
print( 'Everything is done.' )
    
