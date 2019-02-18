## Simple tool to quickly run Karax applications. Generates the HTML
## required to run a Karax app and opens it in a browser.

import os, strutils, parseopt, browsers, strformat, shell

const
  css = """
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bulma/0.5.3/css/bulma.min.css">
  <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css">
"""
  html = """
<!DOCTYPE html>
<html>
<head>
  <title>plotDataClient</title>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <style>
  .dropbtn {
      background-color: #3498DB;
      color: white;
      padding: 16px;
      font-size: 16px;
      border: none;
      cursor: pointer;
  }

  .dropbtn:hover, .dropbtn:focus {
      background-color: #2980B9;
  }

  .dropdown {
      position: relative;
      display: inline-block;
  }

  .dropdown-content {
      display: none;
      position: absolute;
      background-color: #f1f1f1;
      min-width: 160px;
      overflow: auto;
      box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
      z-index: 1;
  }

  .dropdown-content a {
      color: black;
      padding: 12px 16px;
      text-decoration: none;
      display: block;
  }

  .svg-wrap {
      background-color:red;
      height:0;
      padding-top:63.63%; /* 350px/550px */
      position: relative;
  }

  .plot-style{
      width: 60%;
      margin: auto;
  }

  svg {
      height: 100%;
      display:block;
      width: 100%;
  }

  #grid {
      grid-column-gap: 20px;
  }
  #grid {
      display: grid;
      height: 100px;
      grid-template-columns: repeat(10, 1fr);
      grid-template-rows: 100px;
      column-gap: 20px;
  }

  .dropdown a:hover {background-color: #ddd;}

  .show {display: block;}
  </style>
  $2
</head>
<body id="body">
<div id="ROOT" />
<script type="text/javascript" src="$1.js"></script>
</body>
</html>
"""

proc exec(cmd: string) =
  if os.execShellCmd(cmd) != 0:
    quit "External command failed: " & cmd

proc writeDataPath(name: string) =
  ## write the data file to file being included by `staticClient.nim`
  var outfile = open("resources/data_input.nim", fmWrite)
  let content = getAppDir() / name
  outfile.write(&"const fname = \"{content}\"")
  outfile.close()

proc main =
  var op = initOptParser()
  var rest = op.cmdLineRest
  var file = ""
  # set run and css as defaults
  var run = true
  var buildRelease = false
  var selectedCss = css

  while true:
    op.next()
    case op.kind
    of cmdLongOption:
      case op.key
      of "run":
        run = true
        rest = rest.replace("--run ")
      of "css":
        selectedCss = css
        rest = rest.replace("--css ")
      of "file":
        # set the file to be loaded
        rest = rest.replace("--file:")
        rest = rest.replace(op.val)
        writeDataPath(op.val)
      else: discard
    of cmdShortOption:
      case op.key
      of "r":
        run = true
        rest = rest.replace("-r ")
      of "d":
        if op.val == "release":
          buildRelease = true
          rest = rest.replace("-d:")
          rest = rest.replace(op.val)
    of cmdArgument: file = op.key
    of cmdEnd: break

  if file.len == 0: quit "filename expected"
  let name = file.splitFile.name
  createDir("nimcache")
  if buildRelease:
    echo "Building ", name, " in release mode"
    exec("nim js -d:release --out:nimcache/" & name & ".js " & rest)
  else:
    echo "Building ", name, " in debug mode"
    exec("nim js --out:nimcache/" & name & ".js " & rest)
  let dest = "nimcache" / name & ".html"
  writeFile(dest, html % [name, selectedCss])
  if run: openDefaultBrowser(dest)

main()
