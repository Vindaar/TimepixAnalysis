import os, ospaths
import strutils
import times


proc parseTOSDateString*(date_str: string): Time = 
  # function receives a string from a date time from TOS and creates
  # a Time object from it
  let date = toTime(parse(date_str, "yyyy-MM-dd.hh:mm:ss"))
  return date

proc readDateFromEvent*(filepath: string): string = 
  var date_str = ""
  for line in lines filepath:
    if "dateTime" in line:
      # if dateTime is found, split, assign to string and break from while
      let line_seq = split(line, " ")
      date_str = line_seq[high(line_seq)]
      break
    else:
      continue

  return date_str

proc readListOfFilesAndGetTimes*(path: string, list_of_files: seq[string]): seq[Time] = 
  # function analogues to walkRunFolderAndGetTimes except we read all files
  # given by list_of_files
  var i: int = 0
  result = @[]

  for file in list_of_files:
    let filename = joinPath(path, file)
    let date_str = readDateFromEvent(filename)
    if date_str isnot "":
      result.add(parseTOSDateString(date_str))
    if i mod 500 == 0:
      echo i, " files read."
    i += 1    
  return result

proc walkRunFolderAndGetTimes*(folder: string): seq[Time] = 
  var i: int = 0
  result = @[]

  # for file in walkDir(folder):
  #   let path = file.path
  for file in walkFiles(joinPath(folder, "data*.txt-fadc")):
    let path = file
    if "data" in path and "fadc" in path:
      let filename = strip(path, false, true, {'-', 'f', 'a', 'd', 'c'})
      let date_str = readDateFromEvent(filename)
      if date_str isnot "":
        result.add(parseTOSDateString(date_str))
      if i mod 500 == 0:
        echo i, " files read."
      i += 1
  return result

proc writeDateSeqToFile*(date_seq: seq[Time]): void = 
  let outfile = "times.txt"
  var f: File
  if open(f, outfile, fmWrite):
    f.write("# date\t unix time\n")
    for time in date_seq:
      let sec: string = formatFloat(toSeconds(time))
      #let t: string = format(getLocalTime(time), "yyyy-MM-dd hh:mm:ss")
      #let str = join(@[t, sec, "\n"], sep = "\t")
      let str = join(@[$(time), sec, "\n"], sep = "\t")
      f.write(str)
  f.close()

