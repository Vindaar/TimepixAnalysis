#!/usr/bin/env python
# this script converts Christoph's CDL reference ROOT trees
# to HDF5 files

import numpy as np
from ROOT import TTree, TFile
import os
import argparse
import h5py
from helper_functions import get_energy_bins, get_energy_binning, get_hdf5_background_layout, get_hdf5_calibration_layout, get_energy_bin_for_energy

# def get_energy_index(energy, energy_bins):
#     # this function gets the correct index for an event energy
#     # corresponding to the energy binnings
import time

def getall(d, basepath="/"):
    "Generator function to recurse into a ROOT file/dir and yield (path, obj) pairs"
    for key in d.GetListOfKeys():
        kname = key.GetName()
        kclass = key.GetClassName()
        if kclass == "TTree":
            # in case we're dealing with a tree, simply return the tree
            yield kname, d.Get(kname)
        elif key.IsFolder():
            # TODO: -> "yield from" in Py3
            for i in getall(d.Get(kname), basepath+kname+"/"):
                yield i
        else:
            yield basepath+kname, d.Get(kname)

def create_or_get_dataset(col_dict, br_name, data_type, nbins, dtype = None, shape = (10000, 1), maxshape = (None, 1)):
    """
    col_dict : dictionary containing the different trees / files (for the columns, hence the name)
    br_name  : string of the name of the dataset we wish to create or get
    data_type: not to be confused with dtype, data_type is "calibration" or "background" of the data
               we read
    nbins    : (application specific) the number of energy bins we use for our hdf5 file
    dtype    : valid numpy dtype or hdf5 special type (if dtype is None, we do not specify manually)
    shape    : tuple of the shape for the data (similar to numpy shape)
    maxshape : maximum shape the dataset may be enlarged to

    This function first trys to access the dataset with the name 'br_name'. If it does not
    yet exist it is created.
    """
    try:
        dset  = col_dict[data_type][br_name]
    except KeyError:
        # if a KeyError is thrown, the dataset does not exist, thus create it
        if dtype is None:
            print("Now here ", br_name)
            dset  = col_dict[data_type][br_name].create_dataset(br_name,
                                                                shape, maxshape = maxshape,
                                                                dtype=dtype)
    return dset

def determine_ind_counter(dset, nbins):
    """
    dset : a dataset (which has the same size as any of the other dsets
           which one wants to write to), from which ind_counter is det.
    nbins: number of energy bins

    This function determines the ind_counter, meaning the index for each
    dataset, in case the dsets already have data written to them. Necessary,
    since we enlarge dsets by 10000 lines each time, but may only write
    partially (rest still zero)
    """
    # in this case rebuild the ind_counter
    ind_counter = np.zeros(nbins, dtype=np.int)
    for i in xrange(nbins):
        # we create an array of the size of all elements in each
        # bin. x_cols[i][j][0]. The 0 simply gives us the actual array of each
        # event. nonzero returns an array of the non zero elements.
        # nonzero(...)[0][-1] + 1
        # because: [0] to get the actual array
        #          [-1] to get the last element
        #          + 1 to use the next element as our first index
        try:
            ind_counter[i] = np.nonzero([ np.size(dset[i][j][0])
                                          for j in xrange(len(dset[i])) ])[0][-1] + 1
        except IndexError:
            # in this case there is no element yet
            ind_counter[i] = 0

    return ind_counter



def read_from_tree_and_write(tree, data_type):

    # create special hdf5 type
    dt8  = h5py.special_dtype(vlen=np.dtype(np.uint8))
    dt16 = h5py.special_dtype(vlen=np.dtype(np.uint16))


    # TODO: add the following to the output file
    #       - RMS_y = float
    #       - number of pixels = int
    #       - pos_x & pos_y = float & float
    #       - add to a unified dataset with energy? datatypes?


    # create empty buffer lists for the data columns, one element for
    # each list (= one energy bin).
    # once one list (= one energy bin) contains 10.000 events, we write
    # it to file
    nbins = len(energy_binning) - 1
    # x_cols  = create_or_get_dataset("x_col", nbins, dt8)
    # y_cols  = create_or_get_dataset("y_col", nbins, dt8)
    # ch_cols = create_or_get_dataset("ch_col", nbins, dt16)
    # E_cols  = create_or_get_dataset("E_col", nbins, np.float)


    # NOTE: thinkable to do the following to automate this:
    # get all branches of a tree tree.GetListOfBranches()
    # and get all names of these using branch.GetName()
    # create list of names, create dictionary of "name" and a dataset
    # fill automatically by running over all events and

    # TESTING
    print(tree)
    for br in tree.GetListOfBranches():
        print(br.GetName())
    branchNameList = [br.GetName() for br in tree.GetListOfBranches() ]
    # now create datasets for these
    dset_dict = {}
    for branch in branchNameList:
        # going over all branches and create datasets using any datatype, unless
        # we're dealing with x_cols, y_cols and ch_cols
        # we exclude the x, y, charge vectors, since we don't need them for the likelihood
        if branch in ["XCoordinatesVector", "YCoordinatesVector"]:
            # dset_dict[branch] = create_or_get_dataset(col_dict, branch, data_type, nbins, dtype = dt8)
            pass
        elif branch in ["ChargeValuesVector"]:
            # dset_dict[branch] = create_or_get_dataset(col_dict, branch, data_type, nbins, dtype = dt16)
            pass
        else:
            dset_dict[branch] = create_or_get_dataset(col_dict, branch, data_type, nbins)

    # now determine ind_counter (need to hand any element from dictionary
    ind_counter = determine_ind_counter(dset_dict["XCoordinatesVector"], nbins)

    # now that we have the tree, run over all events and split
    print('Total number of entries in this tree: %i' % tree.GetEntries())
    for j, event in enumerate(tree):
        if j % 10000 == 0 and j > 0:
            print('%i events done.' % j)
        # now get all elements of tree for event

        # deal with the following for variables differently, because of special datatypes
        # and Energy needed for binning
        x = y = ch = E = None
        # first get the energy
        E = float(event.EnergyFromCharge)

        # now get correct index for group of energy range
        ind = get_energy_bin_for_energy(E)
        if ind == -1:
            # in this case the energy is smaller than 0.15 keV and thus argmin returns
            # 0. Due to python indexing index -1 represents the last element of a list
            # don't want to add 0 < E < 0.15 to last bin
            # in case of E > 10, argmin also returns 0, thus energies of E > 10, also
            # excluded, as wanted
            continue

        if ind_counter[ind] % 10000 == 0 and ind_counter[ind] > 0:
            oldsize = len(dset_dict["XCoordinatesVector"][ind])
            newsize = oldsize + 10000
            print('resizing energy bin %i from %i to %i' % (ind, oldsize, newsize))
            # resize each element of the dictionary
            for val in dset_dict.itervalues():
                val[ind].resize((newsize, 1))

        # given the bin we have for the energy, add the data as a
        # dataset to the correct groups
        for branch in branchNameList:
            if branch is "XCoordinatesVector":
                x   = np.asarray(event.XCoordinatesVector, dtype=np.uint8)
                dset_dict[branch][ind][ind_counter[ind]]  = x
            elif branch is "YCoordinatesVector":
                y   = np.asarray(event.YCoordinatesVector, dtype=np.uint8)
                dset_dict[branch][ind][ind_counter[ind]]  = y
            elif branch is "ChargeValuesVector":
                ch  = np.asarray(event.ChargeValuesVector, dtype=np.uint16)
                dset_dict[branch][ind][ind_counter[ind]] = ch
            elif branch is "EnergyFromCharge":
                E   = float(event.EnergyFromCharge)
                dset_dict[branch][ind][ind_counter[ind]]  = E
            else:
                # in every other case simply access the correct attribute via its name as a
                # string using a python function
                dset_dict[branch][ind][ind_counter[ind]] = getattr(event, branch)

        # and increase the ind_counter, which is the event number for each
        # energy bin
        ind_counter[ind] += 1


def getData(th1d):
    """
    Returns the content of a TH1D as a numpy array
    """
    return np.asarray([x for x in th1d])


def create_hdf5_file(outfile_str, fname, trees):

    # the way this works:
    # - read from ROOT file (calibration file)
    # - run over each tree (different targets from CDL)
    # - read each dataset, excl vector-like sets and write to H5 with same
    #   name and group structure

    # create hdf5 file
    if outfile_str == None:
        print(fname)
        outfile_str = os.path.basename(fname).replace("root", "h5")
        print(outfile_str)
    outfile = h5py.File('.hdf5', 'w')
    #outtrees = [ TTree("%f_%f" % (energy_binning[i], energy_binning[i+1])) for i in xrange(len(energy_bins))]

    # now use filepaths to get the layout of the groups, we're using in the HDF5 file
    reference_layout = get_hdf5_calibration_layout()

    f = TFile.Open(fname)
    # finally run over all files and trees and add data to groups
    for tree in trees:
        print(f)
        for name, obj in getall(f):
            # given a tree, write the contents of the tree to a h5 file
            read_from_tree_and_write(obj, "reference")


    outfile.close()


def main(args):
    # main function

    parser = argparse.ArgumentParser(description = 'H5 Data plotter')
    parser.add_argument('file',
                        help = "The ROOT file from which to read data")
    parser.add_argument('--outfile',
                        default = None,
                        dest = "outfile",
                        help = "The H5 file to which to write the data.")
    args_dict = vars(parser.parse_args())
    fname = os.path.abspath(args_dict["file"])
    outfile = args_dict["outfile"]

    tree = "/MyXrayAnalysisProcessor"
    trees = [tree]

    # signal_trees_str    = filepaths.get_signal_trees()
    # background_tree_str = filepaths.get_background_tree()
    # trees = [background_tree_str, signal_trees_str]

    # outfile_str = os.path.join(os.path.dirname(filepaths.get_background_path()),
    #                            "energy_binned_data.root")

    create_hdf5_file(outfile, fname, trees)



if __name__=="__main__":
    import sys
    main(sys.argv[1:])
