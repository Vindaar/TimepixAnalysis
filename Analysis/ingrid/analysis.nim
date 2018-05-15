import sequtils
import sugar
import strutils, strformat
import rdstdin
import docopt

const doc = """
Timepix Analysis tool

Usage:
  analysis [options]
  analysis <HDF5file> [options]

Options:
  -h --help              Show this help
  --version              Show version

Documentation:
  This tool provides a unified interface to the Timepix analysis suite. 
  It's a CLI tool, allowing either to call individual analysis procedures
  manually or loading a JSON file containinig commands to be executed. 
  Alternatively one may call one of a set of predefined analysis procedures,
  which combine different parts, e.g.:
  Fe spectrum: 
    - read data from a run
    - write to H5
    - perform geometry calculations on clusters
    - make rough cut on geometric properties and write Fe spectrum to 
      H5 file
  In the future this tool will also provide access to template Python plotting
  scripts to create plots from the data automatically.
"""

proc main() =
  let args = docopt(doc)
  echo args
  
  while true:
    let line = readLineFromStdin("> ")
    echo &"Line read: {line}"
  
  
  
when isMainModule:
  main()
