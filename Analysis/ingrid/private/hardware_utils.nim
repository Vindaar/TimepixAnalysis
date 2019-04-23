#####################################################
#### Procs specifically realted to hardware #########
#####################################################

proc getSeptemHChip*(chipNumber: int): string =
  ## returns the name of a given SeptemH chip
  const names = ["E6 W69",
                 "K6 W69",
                 "H9 W69",
                 "H10 W69",
                 "G10 W69",
                 "D9 W69",
                 "L8 W69"]
  result = names[chipNumber]
