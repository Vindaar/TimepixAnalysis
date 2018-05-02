# this set of functions deals with handling input data

import numpy as np
import cPickle
import h5py

import scipy.misc

from scipy import ndimage
import theano
import theano.tensor as T

import matplotlib.pyplot as plt

# in order to read data from the events taken by TOS, we need to include
# the septemModule
# for that we have included a symbolic link from the septemModule to this 
# folder
from septemModule.septemFiles import create_files_from_path_combined, read_zero_suppressed_data_file, create_filename_from_event_number
import helper_functions


def get_event_type(type_str):
    # this function returns the correct value (0 or 1) depending on
    # a string, which describes the event type
    if type_str in ["calibration", "Calibration", "signal", "Signal"]:
        event_type = 1
    elif type_str in ["background", "Background"]:
        event_type = 0 
    else:
        descr = """%s is not an allowed string to be used for get_event_type().
                   Use calibration, signal or background (w or w/o capital)""" % type_str
        raise NotImplementedError(descr)
        
    return event_type

def event_display_from_non_shared_dataset(start_index, data_x, data_y):
    # print a random pic and check if it matches
    i = 0
    while(i < 50):
        sel_pic = start_index+i
        pic1 = np.reshape(data_x[sel_pic], (128, 128))
        print "pic1 is type: ", data_y[sel_pic]
        #scipy.misc.imsave("event.bmp", pic1)
        imgplot = plt.imshow(data_x[sel_pic].reshape((128,128)))
        imgplot.set_cmap('spectral')
        plt.colorbar()
        plt.show()
        raw_input("waiting...")
        i += 1


def block_mean(ar, fact):
    assert isinstance(fact, int), type(fact)
    sx, sy = ar.shape
    X, Y = np.ogrid[0:sx, 0:sy]
    regions = sy/fact * (X/fact) + Y/fact
    res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
    res.shape = (sx/fact, sy/fact)
    return res

def get_hdf5_dataset_size(group):
    # this function returns the number of entries in the group 'group'
    # of a HDF5 file
    # get dataset for energy (any is good)
    E_col  = np.asarray(group["EnergyFromCharge"]).flatten()
    # number of entries
    print('Info about E_col \n\n\n')
    print np.shape(E_col), np.size(E_col),
    print E_col
    print E_col[0]
    nEntries = np.count_nonzero(E_col)
    return nEntries    


def read_hdf5_and_return_dataset_size(hdf5_dict, data_type):
    # calls get_hdf5_dataset_size internally
    # this function returns the number of elements stored in an 
    # energy bin 'energy' for signal or background data 
    f = h5py.File(hdf5_dict["hdf5_path"], 'r')
    # now get the correct group from which to read the data given the energy
    # given to this function
    group_name = helper_functions.create_group_name_for_energy_bin(data_type, hdf5_dict["energy"])
    # now get the correct group
    group = f[group_name]
    nEntries = get_hdf5_dataset_size(group)
    f.close()
    return nEntries
    

def read_data_from_hdf5(hdf5_dict, data_type, nEvents, downsample_flag = True, pickle_flag = False, start_iter = 0, gold_cut = False):
    # hdf5_dict: contains "hdf5_path" == filename:
    #     filename: name of HDF5 file
    #                     "energy" == energy:
    #     energy: given this energy the correct bin will be chosen to read data from
    # data_type: a string describing 'Background' or 'Calibration' data
    # this function reads data from a HDF5 file and returns a shared dataset

    # NOTE: only rudimentary implementation so far. reading always first elements
    # of the energy bins and no pickling supported yet

    # TODO: implement reading only from gold region
    #define CHIPREGIONS_GOLD_X_MIN 4.5
    #define CHIPREGIONS_GOLD_X_MAX 9.5
    #define CHIPREGIONS_GOLD_Y_MIN 4.5
    #define CHIPREGIONS_GOLD_Y_MAX 9.5
    # for that: read PositionX and PositionY entries from tree
    # and compare position with cut values
    

    # downsampling factor controls how far we downsample
    # e.g. 2 means: 256 * 256 --> 128 * 128
    # note: currently factor needs to divide 256
    if downsample_flag is True:
        downsample_factor = 2
    else:
        downsample_factor = 1

    filename = hdf5_dict["hdf5_path"]
    energy   = hdf5_dict["energy"]

    f = h5py.File(filename, 'r')
    # now get the correct group from which to read the data given the energy
    # given to this function
    group_name = helper_functions.create_group_name_for_energy_bin(data_type, energy)
    # now get the correct group
    group = f[group_name]
    # get dataset for x pixels
    # x_col  = group["x_col"]
    # y_col  = group["y_col"]
    # ch_col = group["ch_col"]
    # E_col  = group["E_col"]
    x_col  = group["XCoordinatesVector"]
    y_col  = group["YCoordinatesVector"]
    ch_col = group["ChargeValuesVector"]
    E_col  = group["EnergyFromCharge"]
    if gold_cut is True:
        x_pos_col = group["PositionX"]
        y_pos_col = group["PositionY"]


    img_size = int(256 / downsample_factor)
    print('Image size: %i x %i' % (img_size, img_size))
    print('getting entry size for group name %s' % group_name)
    nEntries = get_hdf5_dataset_size(group)
    print(nEntries)
    if nEvents > nEntries:
        print("nEvents chosen bigger than tree! Reduced nEvents to %i" % nEntries)
        nEvents = nEntries
    # elif nEvents < buf:
    #     buf = nEvents

    # now that we have our final size for nEvents, create the data numpy array
    data_x   = np.empty((nEvents, img_size*img_size))
    energies = np.empty((nEvents))
    #nFiles = np.ceil(float(nEvents) / float(buf))

    # given the group, we can finally read data
    # starting from 1, because bug caused file to start writing at 1
    # TODO: fix
    #if start_iter == 0:
        # set to 1 in case 0 is selected because of bug described in lines above
        #start_iter = 1

    n_events_read = 0
    ind_hdf5 = start_iter
    ind = 0
    #for i in xrange(start_iter, start_iter + nEvents):
    print 'nasdfas', nEvents
    # TODO: rewrite reading such that we simply use fancy indexing to read
    # whole batch of events in one line, instead of iterating over all events
    # in this case we need to be smart about the potential gold cut. Then
    # need to be applied to numpy array
    # in case events are thrown out, need to refill arrays
    while n_events_read < nEvents:
        # ind is the variable to access the correct elements of the numpy arrays
        # in which the store the *read* data
        # i is the variable we use to access the correct line in the HDF5 dataset
        ind = n_events_read

        if ind % 5000 == 0:
            print("%i events read, current event #%i" % (n_events_read, ind_hdf5))
            # for now read nEvents given i
        x  = np.asarray(x_col[ind_hdf5][0])
        y  = np.asarray(y_col[ind_hdf5][0])
        ch = np.asarray(ch_col[ind_hdf5][0])
        E  = float(E_col[ind_hdf5])

        if gold_cut is True:
            # check if event is within gold region
            xpos = int(x_pos_col[ind_hdf5][0])
            ypos = int(y_pos_col[ind_hdf5][0])
            
            cut = (xpos >= x_gold_min and 
                   xpos <= x_gold_max and 
                   ypos >= y_gold_min and 
                   ypos <= y_gold_max)
            if cut is False:
                # in case the event is not within the gold region, we just
                # skip one event further
                continue

        # now we have all we need to create an event as a numpy array
        event_ar = np.zeros([256, 256])
        event_ar[x, y] = ch

        if downsample_flag is True:
            event_ar = block_mean(event_ar, 2)

    
        # before we set data_x array, we ravel the event_ar array
        event_ar = event_ar.flatten()
        # set data_x array and energy array
        data_x[ind]   = event_ar
        energies[ind] = E

        # increment counter
        n_events_read += 1
        ind_hdf5 += 1


    # we will now return a list containing two lists [ [], [] ] corresponding to
    # [dataset_x, dataset_y]. we can't create shared datasets yet, because we only
    # read signal or background like data, but we would like to have our signal and
    # backkground data shuffled in the test, train and validate samples
    data_y = np.empty(len(data_x), dtype=int)

    # call get_event_type to get the correct event type from the data_type string
    event_type = get_event_type(data_type)
    data_y[:] = int(event_type)
    nEventsRead = np.size(energies)
    print('Read %i events from group: %s' % (nEventsRead, group_name))

    data = [data_x, data_y]
    return data, energies


def read_data_from_tree(filename, treename, nEvents, downsample_flag = True, pickle_flag = False, start_iter = 0):
    # this function reads data from a given ROOT tree
    # and returns a shared dataset
    # downsampling factor controls how far we downsample
    # e.g. 2 means: 256 * 256 --> 128 * 128
    # note: currently factor needs to divide 256
    if downsample_flag is True:
        downsample_factor = 2
    else:
        downsample_factor = 1

    # this list contains all files, which were pickled in this function call.
    # used to return the full dataset (WARNING: might crash if too
    # much data is being read!)
    pickled_list = []

    print filename
    # import TFile from ROOT (done here, since we don't need it elsewhere)
    from ROOT import TFile

    f = TFile.Open(filename)
    tree = f.Get(treename)

    # before we actually read data, let's determine the type of data we will
    # be reading (0 background like, 1 signal like), based on string of treename
    if "calibration" in treename:
        event_type = get_event_type("calibration")
    elif "background" in treename:
        event_type = get_event_type("background")

    img_size = int(256 / downsample_factor)
    buf = 5000
    nEntries = tree.GetEntries()
    if nEvents > nEntries:
        print("nEvents chosen bigger than tree! Reduced nEvents to %i" % nEntries)
        nEvents = nEntries
    elif nEvents < buf:
        buf = nEvents

    # now that we have our final size for nEvents, create the data numpy array
    data_x = np.empty((nEvents, img_size*img_size))
    nFiles = np.ceil(float(nEvents) / float(buf))

    # now consider start_iter. we need to add start_iter to nEvents
    #nEvents += start_iter
    # now we should still get nEvents, but started from start_iter, since we skip all i
    # in the for loop below, if we're below start_iter

    j = 0
    max_charges = []
    for i, event in enumerate(tree):
        if i < start_iter:
            continue
        else:
            x  = event.XCoordinatesVector
            y  = event.YCoordinatesVector
            ch = event.ChargeValuesVector
            x  = np.asarray(x)
            ch = np.asarray(ch)
            max_charges.append(np.max(ch))
            # now we have all we need to create an event as a numpy array
            event_ar = np.zeros([256, 256])
            event_ar[x, y] = ch
            
            if downsample_flag is True:
                event_ar = block_mean(event_ar, 2)
            
            # before we set data_x array, we ravel the event_ar array
            event_ar = event_ar.flatten()
            # set data_x array
            data_x[j] = event_ar
            if j != 0 and (j+1) % buf == 0 or (j+1) == nEvents:
                print("%i number of events read from tree" % (j+1))
                path = 'data/' + treename + '_' + str(j) + '.dat'
                if pickle_flag is True:
                    print("dumping events to file %s" % path)
                    save_file = open(path, 'wb')
                    cPickle.dump(data_x, save_file, -1)
                    save_file.close()
                    pickled_list.append(path)
                    del(data_x)
                    data_x = []
                # else we continue and do not delete array!
                if (j+1) == nEvents:
                    break
            j += 1

    # if we pickle data, we need to read that back to return the full dataset
    if pickle_flag is True:
        for path in pickled_list:
            load_file = open(path, 'r')
            data_x.extend(cPickle.load(load_file))
            load_file.close()

    # now normalize data_x by max_charge:
    #data_x = data_x / max_charge
    #print 'max value data_x ', np.max(data_x)
    #print 'max charge ', np.percentile(max_charges, 70)
    #data_x = data_x / np.percentile(max_charges, 70)
    # print 'max charge now ', np.max(data_x)
    
    # we will now return a list containing two lists [ [], [] ] corresponding to
    # [dataset_x, dataset_y]. we can't create shared datasets yet, because we only
    # read signal or background like data, but we would like to have our signal and
    # backkground data shuffled in the test, train and validate samples
    data_y = np.empty(len(data_x), dtype=int)
    data_y[:] = int(event_type)

    data = [data_x, data_y]
    return data

def read_data_from_run_folder(folder, 
                              nEvents, 
                              event_type,
                              start_iter = 0,
                              downsample_flag = True, 
                              pickle_flag = False):
    # this function reads data from a run folder
    # it calls the same functions used for the event display
    # in the same way as read_data_from_tree it also downsamples and allows cPickled data
    # to be read
    # the argument event_type defines signal or background events. Needs to be handed to
    # the function, since we cannot know from folder name 
    #     background           = 0
    #     calibration / signal = 1

    # TODO: think about implementing start_iter into this function
    # however, for now we want to get the full run
    # >>> not true, we want only the tracking part or rather we need
    # to seperate tracking and non tracking...

    # downsampling factor controls how far we downsample
    # e.g. 2 means: 256 * 256 --> 128 * 128
    # note: currently factor needs to divide 256
    if downsample_flag is True:
        downsample_factor = 2
    else:
        downsample_factor = 1
    # use factor to determine final image size
    img_size = int(256 / downsample_factor)

    # the first thing to do, before we can read the data, is to create a set of 
    # events, which exist in the run folder
    eventSet, fadcSet = create_files_from_path_combined(folder,
                                                        eventSet,
                                                        fadcSet,
                                                        False)
    # decrease nEvents to length of eventSet if larger
    if nEvents > len(eventSet):
        nEvents = len(eventSet)
    # now create a numpy array of frames for nEvents
    data_x = np.empty((nEvents, img_size*img_size))
    
    for i, eventNumber in enumerate(eventSet):
        if eventNumber < start_iter:
            # if event number smaller than start_iter, continue 
            continue
        elif eventNumber + start_iter > nEvents:
            # in case we have the number of events, we wish to read
            # break from the loop
            break
        else:
            # first create the filename
            filename     = create_filename_from_event_number(eventSet,
                                                             eventNumber,
                                                             nfiles,
                                                             fadcFlag = False)
            # we create the event and chip header object
            evHeader, chpHeaderList = read_zero_suppressed_data_file(folder + filename)
            # now go through each chip header and add data of frame
            # to chip_arrays

            # TODO: need to implement differentiation between tracking and non-tracking datasets

            for chpHeader in chpHeaderList:
                # get data of the frame for this chip
                chip_data = chpHeader.pixData
                # TODO: implement what to do if more than one chip used. Will be necessary, once
                # we need to analyze septem board events
                
                # now create temporary event array, downsample and assign to data_x
                event_ar = np.zeros([256, 256])
                event_ar[chip_data[:,0], chip_data[:,1]] = chip_data[:,2]
                
                if downsample_flag is True:
                    event_ar = block_mean(event_ar, downsample_factor)
                
                # before we set data_x array, we ravel the event_ar array
                event_ar = event_ar.flatten()
                # set data_x array
                data_x[i] = event_ar

            # some output
            if i % 1000 == 0:
                print event_num, ' events done.' 

    # we will now return a list containing two lists [ [], [] ] corresponding to
    # [dataset_x, dataset_y]. we can't create shared datasets yet, because we only
    # read signal or background like data, but we would like to have our signal and
    # backkground data shuffled in the test, train and validate samples
    data_y = np.empty(len(data_x), dtype=int)
    data_y[:] = int(event_type)

    data = [data_x, data_y]
    return data


def shared_dataset(data_xy, borrow=True):
    """ Function that loads the dataset into shared variables
    The reason we store our dataset in shared variables is to allow
    Theano to copy it into the GPU memory (when code is run on GPU).
    Since copying data into the GPU is slow, copying a minibatch everytime
    is needed (the default behaviour if the data is not in a shared
    variable) would lead to a large decrease in performance.
    """
    data_x, data_y = data_xy
    shared_x = theano.shared(np.asarray(data_x,
                                           dtype=theano.config.floatX),
                             borrow=borrow)
    # shared_y = theano.shared(np.asarray(data_y,
    #                                        dtype=theano.config.floatX),
    #                          borrow=borrow)
    shared_y = theano.shared(np.asarray(data_y,
                                           dtype='int32'),
                             borrow=borrow)
    # When storing data on the GPU it has to be stored as floats
    # therefore we will store the labels as ``floatX`` as well
    # (``shared_y`` does exactly that). But during our computations
    # we need them as ints (we use labels as index, and if they are
    # floats it doesn't make sense) therefore instead of returning
    # ``shared_y`` we will have to cast it to int. This little hack
    # lets ous get around this issue
    return shared_x, shared_y#T.cast(shared_y, 'int32')

    # change how?!?!

def get_shared_datasets(file_tree_dict, hdf5_dict, nEvents, pickle_flag, downsample_flag):
    # this function reads data either from a ROOT tree or a HDF5 file
    # either file_tree_dict needs to be given as None or hdf5_dict as None

    if file_tree_dict is not None:
        # now we load some background data
        background_data = read_data_from_tree(file_tree_dict["f_background"],
                                              file_tree_dict["t_background"],
                                              nEvents, 
                                              downsample_flag = downsample_flag,
                                              pickle_flag = pickle_flag)
        # and some signal data
        signal_data     = read_data_from_tree(file_tree_dict["f_signal"],
                                              file_tree_dict["t_signal"],
                                              nEvents, 
                                              downsample_flag = downsample_flag,
                                              pickle_flag = pickle_flag)
    elif hdf5_dict is not None:
        # NOTE: only rudimentary implementation so far. reading always first elements
        # of the energy bins and no pickling supported yet

        # now load some background data
        background_data, b_energies = read_data_from_hdf5(hdf5_dict, 
                                                          "Background",
                                                          nEvents,
                                                          downsample_flag = downsample_flag)
        # and some signal data
        signal_data, s_energies = read_data_from_hdf5(hdf5_dict, 
                                                      "Calibration",
                                                      nEvents,
                                                      downsample_flag = downsample_flag)
    else:
        raise NotImplementedError

    # get seperate parts of both sets
    background_data_x, background_data_y = background_data
    signal_data_x, signal_data_y         = signal_data

    # create set of shuffled indices
    nInd         = np.size(background_data_x[:,0])
    nPix         = np.size(background_data_x[0])
    shuffled_ind = np.random.permutation(2*nInd)
    # now create an empty array, which will store both signal and background
    data_x = np.empty((2*nInd, nPix))
    data_y = np.empty((2*nInd))
    # and fill it
    data_x[:nInd] = background_data_x
    data_x[nInd:] = signal_data_x
    data_y[:nInd] = background_data_y
    data_y[nInd:] = signal_data_y

    # now shuffle the array
    data_x = data_x[shuffled_ind]
    data_y = data_y[shuffled_ind]
    
    # now seperate data into test, valid and evaluate samples
    n_train_index = nInd
    n_test_and_valid_ind = n_train_index + int(nInd / 2.0)
    train_x = data_x[:n_train_index]
    train_y = data_y[:n_train_index]
    valid_x = data_x[n_train_index:n_test_and_valid_ind]
    valid_y = data_y[n_train_index:n_test_and_valid_ind]
    test_x  = data_x[n_test_and_valid_ind:]     
    test_y  = data_y[n_test_and_valid_ind:] 
    # create sets and shared
    train_set = shared_dataset([train_x, train_y])
    valid_set = shared_dataset([valid_x, valid_y])
    test_set  = shared_dataset([test_x,  test_y])

    
    return [train_set, valid_set, test_set]


def get_evaluation_set(hdf5_dict, 
                       root_dict,
                       data_type,
                       nEvents,
                       downsample_flag = True, 
                       pickle_flag = False, 
                       start_iter = 0,
                       nFitInMemory = 30000):
    # this function returns a shared dataset of type eval_type
    # containing nEvents events
    # data_type: string describing the main data type we read from (in case
    #            we only wish to read either signal or background
    #            if None we read both
    # if filename2 and treename2 are set, this function allows to read a second set
    # of data (as non shared variable), which will then be added to the notFit array
    
    if hdf5_dict is not None:
        if nEvents > nFitInMemory:
            descr = """For HDF5 files we do not support reading more than what fits into
                       GPUs memory. Instead we read the next batch to be put into GPU memory
                       on the fly in the evaluation loop. Simplifies code and should not be 
                       much slower."""
            raise NotImplementedError(descr)

        if data_type is None:
            # in case we read both kinds of data, half number of reading for each
            nEvents /= 2
            data, energies   = read_data_from_hdf5(hdf5_dict, 
                                                   "signal", 
                                                   nEvents, 
                                                   downsample_flag, 
                                                   False, 
                                                   start_iter = start_iter)
            data2, energies2 = read_data_from_hdf5(hdf5_dict, 
                                                   "background", 
                                                   nEvents, 
                                                   downsample_flag, 
                                                   False, 
                                                   start_iter = start_iter)
        else:
            # set 
            data2 = None
            data, energies   = read_data_from_hdf5(hdf5_dict, 
                                                   data_type, 
                                                   nEvents, 
                                                   downsample_flag, 
                                                   False, 
                                                   start_iter = start_iter)
    elif root_dict is not None:
        filename = root_dict["filename"]
        treename = root_dict["treename"]
        if len(root_dict) > 2:
            filename2 = root_dict["filename2"]
            treename2 = root_dict["treename2"]
        else:
            filename2 = ''
            treename2 = ''

        # read the correct data from file:
        print('Starting to read data...')
        data = read_data_from_tree(filename, 
                                   treename, 
                                   nEvents, 
                                   downsample_flag = downsample_flag,
                                   pickle_flag = pickle_flag, 
                                   start_iter = start_iter)
        
        data2 = None
        if filename2 is not '' and treename2 is not '':
            data2 = read_data_from_tree(filename2, 
                                        treename2, 
                                        nEvents, 
                                        downsample_flag = downsample_flag,
                                        pickle_flag = pickle_flag, 
                                        start_iter = start_iter)

    # now we need to be careful. If we have more than roughly
    # 30000 events, it won't fit on the GPU anymore. So we need
    # to create only the very first dataset as shared. We will
    # return the rest as a not yet shared dataset. That will be
    # handled during evaluation (change shared datasets in place on
    # GPU)
    if len(data[0]) > nFitInMemory:
        dataFit    = [ data[0][:nFitInMemory], data[1][:nFitInMemory] ]
        dataNotFit = [ data[0][nFitInMemory:], data[1][nFitInMemory:] ]
        if data2 is not None:
            print np.shape(dataNotFit[0]), np.shape(data2[0])
            dataNotFitX = np.concatenate( (dataNotFit[0], data2[0]), axis=0)
            dataNotFitY = np.concatenate( (dataNotFit[1], data2[1]), axis=0)
            dataNotFit = [dataNotFitX, dataNotFitY]
            # slightly complicated to deal with not fitting arrays, if we need to extend them
            nEvents = np.shape(dataFit[0])[0] + np.shape(dataNotFit[0])[0]
            print 'nevents', nEvents


            # len1 = len(dataNotFit[0])
            # len2 = len(data2[0])
            # dataNotFitX = np.empty( ((len1 + len2), np.shape(dataNotFit[0])[1] ) )
            # dataNotFitX = np.empty( (len1 + len2) )
            # dataNotFitX[:len1] = dataNotFit[0]
            # print np.shape(dataNotFit[0]), np.shape(dataNotFitX) 
            # dataNotFitY[:len1] = dataNotFit[1]
            # print np.shape(dataNotFit[1]), np.shape(dataNotFitY)
            # #import sys
            # #sys.exit()
            # dataNotFitX[len1:] = data2[0]
            # dataNotFitY[len1:] = data2[1]
            # dataNotFitX = np.asarray(dataNotFitX)
            # dataNotFitX = np.asarray(dataNotFitY)
            # dataNotFit = [ dataNotFitX, dataNotFitY ]
            #nEvents *= 2
        # and create the shared dataset
        dataset = shared_dataset(dataFit)
        print('...finished reading data, returning only partly as shared dataset')
        return nEvents, dataset, dataNotFit
    else:
        # if we have data2, combine both datasets
        if data2 is not None:
            data_x = np.concatenate( (data[0], data2[0]), axis=0)
            data_y = np.concatenate( (data[1], data2[1]), axis=0)
            # concatenate energy arrays as well
            energies = np.concatenate( (energies, energies2), axis=0 )

            data = [data_x, data_y]
            

        # and create the shared dataset
        dataset = shared_dataset(data)

        print ('...finished reading data')
        if data2 is not None:
            # multiply nEvents in this case by 2, since we divided by 2 before
            nEvents *= 2
        #nEvents = np.size(data[0])
        return nEvents, dataset, energies
            

def get_evaluation_set_from_run_folder(folder, 
                                       nEvents,
                                       downsample_flag = True, 
                                       pickle_flag = False, 
                                       start_iter = 0,
                                       nFitInMemory = 30000):
    # this function creates a shared dataset from files, which are
    # part of a run folder
    # it calls the read_data_from_run_folder function and packs it into
    # a shared dataset afterwards

    # TODO: write wrapper functions around a general 'get_evaluation_set' function
    # which only differentiate what kind of data is read.
    print('Starting to read data...')

    data = read_data_from_run_folder(folder, 
                                     nEvents, 
                                     start_iter,
                                     downsample_flag, 
                                     pickle_flag)
    
    # and create the shared dataset
    dataset = shared_dataset(data)
    
    print ('... finished reading data')
    return nEvents, dataset
