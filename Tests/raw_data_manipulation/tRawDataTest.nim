#import ingrid / ingrid_types
#import ingrid / tos_helpers
#import ingridDatabase / databaseDefinitions
import sequtils, strutils, os, algorithm, strformat
import unittest
#import json, sugar
#from ingrid / reconstruction import recoEvent

import ggplotnim
import seqmath

const pwd = currentSourcePath().parentDir
const dataPwd = pwd / "../../resources/TPAresources/raw_data_manipulation/"
const jsonData = dataPwd / "marlinEvents.json"

import times
let plotSuffix = $getTime().toUnix & ".pdf"

suite "raw data manipulation":
  ## these tests check whether the raw data manipulation produces HDF5
  ## files as we expect them given certain command line arguments
  # first run raw data manipulation
  shell:
    "../../Analysis/raw_data_manipulation" `
