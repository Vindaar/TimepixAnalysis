import os, ospaths
import strutils
import times
import algorithm
import re
import tables

type 
  # an object, which stores information about a run's start, end and length
  RunTimeInfo* = object
    t_start*: Time
    t_end*: Time
    t_length*: TimeInterval


proc parseTOSDateString*(date_str: string): Time = 
  # function receives a string from a date time from TOS and creates
  # a Time object from it
  let date = toTime(parse(date_str, "yyyy-MM-dd.hh:mm:ss"))
  return date

proc readDateFromEvent*(filepath: string): string = 
  # procedure opens the given file and returns the dateTime value
  # (TOS syntaxed date string)
  for line in lines filepath:
    if "dateTime" in line:
      # if dateTime is found, split, assign to string and break from while
      let line_seq = split(line, " ")
      result = line_seq[high(line_seq)]
      break
    else:
      continue

proc formatAsOrgDate*(t: Time, org_format = "yyyy-MM-dd ddd H:mm"): string =
  # this procedure formats the given Time object as an org-mode date
  # by first converting it to a TimeInfo object
  let ti = getLocalTime(t)
  result = format(ti, org_format)

proc getTimeFromEvent*(file: string): Time = 
  # this procedure returns the time info from an event,
  # given the filename by parsing the TOS date syntax
  let date_str = readDateFromEvent(file)
  result = parseTOSDateString(date_str)


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


proc getRunTimeInfo*(run_files: seq[string]): RunTimeInfo =
  # this procdure creates a RunTimeInfo object from a given list of files
  # in a run folder (not sorted). The list is sorted and then the first and
  # last event are read, converted to Time objects and the length of the run 
  # is calculated from the difference between both events
  result = RunTimeInfo()

  # sort list of files
  let sorted_files = sorted(run_files, system.cmp[string])
  let first_file = sorted_files[0]
  let last_file  = sorted_files[^1]

  let time_first = getTimeFromEvent(first_file)
  let time_last  = getTimeFromEvent(last_file)

  let run_length = initInterval(seconds=int(time_last - time_first))
  
  result.t_start = time_first
  result.t_end = time_last 
  result.t_length = run_length


# type
#   EventHeader* = object
    


# proc readEventHeader*(filepath: string): Table[string, string] =
#   # this procedure reads a whole event header and returns 
#   # a table containing the data, where the key is the key from the data file

#   # we define a regex for the header of the file
#   let regex = r"^##.(.*:)\S(.*)$"
#   for line in lines filepath:
#     echo line, "\t", line.match(re(regex))
    #if line.match(re):
    #   # if dateTime is found, split, assign to string and break from while
    #   let line_seq = split(line, " ")
    #   result = line_seq[high(line_seq)]
    #   break
    # else:
    #   continue

  
