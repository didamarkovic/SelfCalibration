# Rory Holmes
# Mar 2012

# This file contains a sample parameter set for the self-calibration
# simulations. All the parameters are stored in a single dictionary.
# This file is parsed with the eval() function.


{
    # Programmatic Parameters
    # =======================
    # String: The output directory path for the simulation run
    'data_dir'            :     'default_output', 
    # Boolean: Set to True to save out all the data required for plots
    'plotdata'            :     True,
    # Boolean: Set to True to run the simulation in verbose mode,
    'verbose'             :     False, 

    # Sky Parameters
    # ==============
    # Float: the total number of sources (all magnitude) per unit area
    # on the sky
    'density_of_stars'    :     50., 
    # Float array: The parameters describing the magnitude distribution
    # of the sources in the sky, according to
    # log10(dN/dm) = A + B * mag + C * mag ** 2
    'powerlaw_constants'  :     [-13.34863146, 1.25429311, -0.02122949], 
    # Float Array: The area of sky to generate sources in 
    # [alpha_min, alpha_max, beta_min, beta_max]
    'sky_limits'          :     [-4.0, 4.0, -4.0, 4.0],

    # Instrument Parameters
    # =====================
    # Float Array: The simulate imager's field-of-view [alpha, beta]
    'FoV'                 :     [0.76, 0.72],
    # Float: The saturation limit of the simulated imager
    'm_min'               :     17.,
    # Float: The 10-sigma detection limit of the simulated imager
    'm_max'               :     22.,   
    # Floats: The parameters used in the measurement noise model
    # sigma ** 2 = (1 + epsilon) * delta ** 2 + eta ** 2 * count_rate ** 2
    'eta'                 :     0.00173214,
    'delta'               :     0.1584655417,
    # where epsilon is a random number drawn uniformly in the range
    # [0.0, epsilon_max) for each measurement
    'epsilon_max'         :     1.,

    # Simulation Parameters
    # =====================
    # String: The file path to the survey strategy to be performed
    # during the simulation run.
    'survey_strategies'   :     ['survey'],
    # Float: The fraction of the bright sources to be used by the
    # self-calibration simulations. Reduce to increase simulation time.
    'useful_fraction'     :     1.,
    # Integer Array: The number of sample points on the focal plane for the
    # best-in-basis fitting and all of the badness measures
    'ff_samples'          :     [300, 300],
    # Integer: The order of the flat-field used to fit the instrument response
    # in the self-calibration procedure
    'flat_field_order'    :     8,
    # Float: The seed for the random number generator (to ensure repeat runs
    # get the same answer!)
    'seed'                :     1.,
    # Float: The stop condition for the self-calibration procedure and the
    # best-in-basis fitting (stop when difference is less than 2 times this)
    'stop_condition'	    :     1e-5,
    # Integer: The maximum number of iterations in the self-calibration
    # procedure and the best-in-basis fitting
    'max_iterations'	    :     1049
    } 
