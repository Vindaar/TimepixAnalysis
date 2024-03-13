import parseutils, strutils, os, sugar
import cligen


proc main(input: string, outfile = "cdl_cuts_") =
  let outfile = if outfile == "cdl_cuts_": outfile & input.extractFilename.split('-')[^1] else: outfile
  var
    tf: string
    muStr: string
    sStr: string
    nStr: string
    mu: float
    sigma: float
    N: float
    tfSt = "svg:"
    muSt = "gmu = "
    sSt = "gs  = "
    nSt = "gN  = "
  var outf = open(outfile, fmWrite)
  for line in lines(input):
    if line.startsWith(tfSt):
      tf = line.dup(removePrefix(tfSt))
    elif line.startsWith(muSt):
      discard parseUntil(line, muStr, until = {'\\'}, start = muSt.len)
      mu = muStr.strip.parseFloat
    elif line.startsWith(nSt):
      discard parseUntil(line, nStr, until = {'\\'}, start = muSt.len)
      N = nStr.strip.parseFloat
    elif line.startsWith(sSt):
      discard parseUntil(line, sStr, until = {'\\'}, start = sSt.len)
      sigma = sStr.strip.parseFloat
      # write to file
      outf.write(tf & ": " & $(mu - 3 * sigma) & ", " & $(mu + 3 * sigma) & " at amplitude: " & $N & "\n")
  outf.close()

when isMainModule:
  dispatch main
