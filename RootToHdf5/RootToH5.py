#!/usr/bin/env python

# A soon to be general purpose ROOT to HDF5 converter
# Currently still a lot of hardcoded stuff and also it's extremely slow

import numpy as np
from ROOT import TTree, TFile
import os
import argparse
import h5py

import time

def getall(d, basepath="/"):
    "Generator function to recurse into a ROOT file/dir and yield (path, obj) pairs"
    for key in d.GetListOfKeys():
        kname = key.GetName()
        kclass = key.GetClassName()
        print(kname)
        if kclass == "TTree":
            # in case we're dealing with a tree, simply return the tree
            yield kclass, kname, key.ReadObj()#d.Get(kname)
        elif key.IsFolder():
            # TODO: -> "yield from" in Py3
            for i in getall(key.ReadObj(), basepath+kname+"/"):#getall(d.Get(kname), basepath+kname+"/"):
                yield i
        else:
            yield kclass, basepath+kname, key.ReadObj()#d.Get(kname)

def getData(th1d, bins = None):
    """
    Returns the content of a TH1D as a numpy array
    """
    if bins is None:
        return np.asarray([x for x in th1d])
    else:
        return np.asarray([th1d.GetBinContent(b) for b in xrange(len(bins))])

def read_from_dir_and_write(h5file, obj, data_type, kclass, fname):
    """
    Equivalent to the read_from_tree_and_write proc except working on
    dictionaries containing TH1Ds instead of leavess of a TTree
    """
    print(kclass, obj)
    if kclass == "TH1D":
        print(kclass)
        nbins = obj.GetNbinsX()
        bins = np.asarray([obj.GetBinCenter(i) for i in xrange(nbins)])
        # get data for each bin
        data = getData(obj, bins)
        # stack them to two columns (x, y)
        data = np.column_stack([bins, data])
        # modify the name of the dataset in case it contains '/'
        kname = obj.GetName()
        br_name = os.path.join(fname, kname)
        dset = h5file.create_dataset(br_name, np.shape(data))
        dset[:] = data


def read_from_tree_and_write(h5file, tree, data_type):

    # TODO: add the following to the output file
    #       - RMS_y = float
    #       - number of pixels = int
    #       - pos_x & pos_y = float & float
    #       - add to a unified dataset with energy? datatypes?

    branchSizeDict = {}
    for br in tree.GetListOfBranches():
        branchSizeDict[br.GetName()] = br.GetEntryNumber()
    print(branchSizeDict)
    # create special hdf5 type
    dt8  = h5py.special_dtype(vlen=np.dtype(np.uint8))
    dt16 = h5py.special_dtype(vlen=np.dtype(np.uint16))

    # now create datasets for these
    dset_dict = {}
    for br in tree.GetListOfBranches():
        name = br.GetName()
        br_name = os.path.join(tree.GetName(), name)
        shape = (branchSizeDict[name], 1)
        if name in ["XCoordinatesVector", "YCoordinatesVector"]:
            dset_dict[name] = h5file.create_dataset(br_name, shape, dtype = dt8)
        elif name in ["ChargeValuesVector"]:
            dset_dict[name] = h5file.create_dataset(br_name, shape, dtype = dt16)
        else:
            dset_dict[name] = h5file.create_dataset(br_name, shape)

    for ind, event in enumerate(tree):
        if ind % 10000 == 0:
            print(ind, " events done")
        for br in tree.GetListOfBranches():
            branch = br.GetName()
            if branch == "XCoordinatesVector":
                x   = np.asarray(event.XCoordinatesVector, dtype=np.uint8)
                dset_dict[branch][ind]  = x
            elif branch == "YCoordinatesVector":
                y   = np.asarray(event.YCoordinatesVector, dtype=np.uint8)
                dset_dict[branch][ind]  = y
            elif branch == "ChargeValuesVector":
                ch  = np.asarray(event.ChargeValuesVector, dtype=np.uint16)
                dset_dict[branch][ind] = ch
            elif branch == "EnergyFromCharge":
                E   = float(event.EnergyFromCharge)
                dset_dict[branch][ind]  = E
            else:
                # in every other case simply access the correct attribute via its name as a
                # string using a python function
                dset_dict[branch][ind] = getattr(event, branch)


def create_hdf5_file(outfile_str, fname):
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
    outfile = h5py.File(outfile_str, 'a')

    f = TFile.Open(fname)
    for kclass, name, obj in getall(f):
        # given a tree, write the contents of the tree to a h5 file
        if kclass == "TTree":
            read_from_tree_and_write(outfile, obj, "reference")
        else:
            read_from_dir_and_write(outfile, obj, "reference", kclass, fname)

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

    create_hdf5_file(outfile, fname)

if __name__=="__main__":
    import sys
    main(sys.argv[1:])
