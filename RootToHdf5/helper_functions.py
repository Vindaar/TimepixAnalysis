# this file contains several functions, which serve as helper functions
import numpy as np
import os
import handle_data

def get_current_best_models():
    # this function contains a list of the current best models for 
    # each energy range (the filename of the trained model)
    # NOTE: the model for 1.2 is still using 128 x 128 images!!

    models = ["unfinishedNetworks/classifier_experimental_cnn_hdf5_0.3_best_test_12.16.cnn",
              #"unfinishedNetworks/classifier_experimental_cnn_hdf5_0.6_best_test_13.98.cnn",
              "unfinishedNetworks/classifier_experimental_ReLU_bigger_cnn_hdf5_0.6_best_test_12.3333333333.cnn",
              #"unfinishedNetworks/classifier_experimental_cnn_hdf5_1.2_best_test_9.8.cnn",
              "unfinishedNetworks/classifier_experimental_cnn_hdf5_1.2_best_test_8.25.cnn",
              "unfinishedNetworks/classifier_experimental_cnn_hdf5_1.5_best_test_3.62.cnn",
              "unfinishedNetworks/classifier_experimental_cnn_hdf5_2.5_best_test_3.58.cnn",
              #"unfinishedNetworks/classifier_experimental_cnn_hdf5_4.0_best_test_7.34.cnn",
              "unfinishedNetworks/classifier_experimental_cnn_hdf5_4.1_best_test_3.33333333333.cnn",
              "unfinishedNetworks/classifier_experimental_cnn_hdf5_5.0_best_test_3.46.cnn",
              "unfinishedNetworks/classifier_experimental_cnn_hdf5_8.0_best_test_5.24.cnn"]
    return models

def get_energy_binning():
    # returns the energy binning boundaries
    energy_binning = np.asarray([0.15,
                                 0.4,
                                 0.7,
                                 1.2,
                                 2.1,
                                 3.2,
                                 4.9,
                                 6.9,
                                 10.0])
    return energy_binning

def get_energy_bins():
    # returns a list of tuples, which give the actual bins
    # of the energy binning
    energy_bins = [ (0.15, 0.4), 
                    (0.4, 0.7), 
                    (0.7, 1.2), 
                    (1.2, 2.1),
                    (2.1, 3.2),
                    (3.2, 4.9),
                    (4.9, 6.9),
                    (6.9, 10.0) ]
    return energy_bins

def get_energy_bin_for_energy(E):
    # this function returns the correct bin index for a given energy E
    # corresponding to the energy bins defined in get_energy_bins()
    # given np.argmin gives us the first element bigger than E,
    # we choose ind - 1 as our index for the correct group

    energy_binning = get_energy_binning()
    ind = np.argmin(energy_binning < E) - 1
    return ind

def get_energy_from_name(name):
    # this function returns the energy a classifier was trained on from 
    # the name of the classifier (since it's included in the name as
    #     "experimental_cnn_hdf5_%.1f" % hdf5_dict["energy"]
    # NOTE: deprecated. from now on we add hdf5_dict to info_dict of classifier
    energy = name.split("_")[-1]
    if '.txt' in energy:
        energy = energy.split('.txt')[0]
    return float(energy)

def get_energy_tuple_from_name(name):
    # this function returns the tuple of the energy bin boundaries
    # based on the name of the classifier
    energy_bins = get_energy_bins()
    energy = get_energy_from_name(name)
    e_ind = get_energy_bin_for_energy(energy)
    e_bin = energy_bins[e_ind]

    return e_bin

def get_hdf5_calibration_layout():
    return "Events/Calibration/E_%.2f_%.2f/"

def get_hdf5_background_layout():
    return "Events/Background/E_%.2f_%.2f/"
    
def create_group_name_for_energy_bin(data_type, E):
    # this function creates the correct group name for 
    # a given energy corresponding to a energy bin in the
    # hdf5 file
    bins = get_energy_bins()
    ind  = get_energy_bin_for_energy(E)
    en_tuple = bins[ind]
    
    if data_type in ["calibration", "Calibration", "signal", "Signal"]:
        layout = get_hdf5_calibration_layout()
    elif data_type in ["background", "Background"]:
        layout = get_hdf5_background_layout()
    else:
        raise NotImplementedError('%s is not a valid data type' % data_type)
        
    name = layout % (en_tuple[0], en_tuple[1])
    return name

def get_classifier_prefix():
    return "p_y_given_x"

def get_classifier_folder():
    return "classifier_eval_outputs"

def get_sig_vs_back_folder():
    return "p_y_given_x_sig_vs_back"

def get_evaluate_output_filename(classifier_fname, sig_vs_back = False):
    # given a filename of a stored classifier, create the 
    # correct filename and path to store the evaluation output
    
    basename  = os.path.basename(classifier_fname)
    # TODO: get rid of the following 'hack'
    if 'best' in basename:
        basename = basename.split("_best")[0].split('classifier_')[-1]

    new_fname = get_classifier_prefix() + "_" + basename + ".txt"

    if sig_vs_back == True:
        # in case of using signal + background evaluation data, add correct folder
        new_fname = os.path.join(get_sig_vs_back_folder(), new_fname)

    classifier_filename = os.path.join(get_classifier_folder(), new_fname)

    return classifier_filename

def read_classifier_evaluation_file(cl_fname):
    # read the data from a classifier p_y given x output file
    outputs     = []
    exp_outputs = []
    
    with open(cl_fname) as infile:
        for line in infile:
            if "#" not in line:
                # ignore the header of the file
                line = line.split()
                outputs.append( float(line[0]) )
                exp_outputs.append( float(line[-2]) )

    outputs     = np.asarray(outputs, dtype = np.float)
    exp_outputs = np.asarray(exp_outputs, dtype = np.float)
    
    # zip the output and return
    out = zip(outputs, exp_outputs)
    
    return out

def calc_efficiency(p_y_given_x, data_type, cut_vals):
    # p_y_given_x is an output from read_classifier_evaluation_files
    # a list of zipped tuples from (output, exp_output)

    # get correct data_type for output (whether 0 or 1 represents a signal)
    event_type = handle_data.get_event_type(data_type)

    # start by getting the correct elements from list (either signal or background)
    outs = []

    for tup in p_y_given_x:
        output, exp_output = tup
        if exp_output == event_type:
            outs.append(output)

    # given a list of signal outputs, calculate signal efficiency 
    # based on integral in range 0 to 1
    # in order to accelerate calculations, we sort signal array
    outs = sorted(outs)

    def efficiency(ar, val, n_tot):
        # calculate efficiency for value val of array ar by finding index 
        # of first value bigger than val
        n     = float(np.size(np.where(ar < val)[0]))
        eff   = n / n_tot
        return eff

    # get number of classified elements
    n_tot = float(np.size(outs))
    effs = np.asarray([efficiency(outs, cut, n_tot) for cut in cut_vals])
    
    background_type = handle_data.get_event_type("background")
    if event_type == background_type:
        # in case of background we're interested in the rejection, not eff,
        # thus subtract from 1
        effs = 1 - effs
    
    return effs


def calc_cut_value(model, data_type, eff):
    # this function calculates the cut value to achieve 
    # an efficiency of eff for a given data type 
    # (signal efficiency for "signal" or background rej. for 
    # "background") for a given model

    # # create classifier output filenames
    fname = get_evaluate_output_filename(model, True)
    
    # read corresponding file
    p_y_given_x = read_classifier_evaluation_file(fname)

    # define cut values 
    cut_vals = np.linspace(0, 1, 1000)

    # use calc_efficiency to calculate the efficiencies
    eff_per_cut = calc_efficiency(p_y_given_x, data_type, cut_vals)

    # given eff_per_cut, find the cut index which corresponds to eff
    cut_ind = np.argmax(eff_per_cut > eff)
    # and use this index and cut_vals to determine cut_val
    cut_val = cut_vals[cut_ind]

    return cut_val


def get_events_passing_cut(model, cut_val):
    # this function returns a list of all events (their energies to be specific)
    # of events which pass the cut_val for a given model
    
    # given model, get the correct filename for the file containing all classified events
    print('Model %s' % model)
    fname = get_evaluate_output_filename(model, sig_vs_back = False)

    events = []
    with open(fname, 'r') as infile:
        for line in infile:
            if "#" not in line:
                line = line.split()
                P_back = float(line[0])
                energy = float(line[-1])
                if P_back < cut_val:
                    events.append(energy)

    return events
                
