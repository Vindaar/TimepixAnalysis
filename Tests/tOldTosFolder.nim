import ingrid/tos_helpers
import os

proc main() =

  const path = "/mnt/Daten/CAST/2014_15/DataRuns/100-Run140613_16-47-30"

  let (a, b, c, d) = isTosRunFolder(path)
  echo "is rf ", a
  echo "runNumber ", b
  echo "run folder kind ", c
  echo "contains rf ", d

  doAssert b == 100 

when isMainModule:
  main()
