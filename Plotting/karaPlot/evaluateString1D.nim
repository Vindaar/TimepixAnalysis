import sequtils, seqmath, strformat, os, shell, strutils, json

let a = 5.5
let b = 12.5

let low = 0.0
let high = 150.0

let fun = &"{a} * sin({a} * x) + {b}"

let tmpl = """
import sequtils, json, seqmath, math
let low = $1
let high = $2
func test(x: float): float =
  result = $3
let res = linspace($1, $2, 30).mapIt(test(it))

let jsonOut = %* {
  "func" : "$3",
  "low" : low,
  "high" : high,
  "res" : res
}
echo jsonOut.pretty
""" % [$low, $high, $fun]

var f = open("/tmp/newFunc.nim", fmWrite)
f.write(tmpl)
f.close()

setCurrentDir(getAppDir())
var funRes = ""
shellAssign:
  funRes = nim "c --hints:off -r /tmp/newFunc.nim"

var output = funRes.parseJson

echo output.pretty
