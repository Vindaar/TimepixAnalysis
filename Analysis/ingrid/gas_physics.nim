import pkg / [xrayAttenuation, unchained, datamancer]

proc absorptionLengthCAST*(energy: keV): cm =
  ## Returns the absorption length for the given energy in keV for CAST
  ## gas conditions:
  ## - Argon / Isobutane 97.7 / 2.3 % ( NOTE : currently Isobutane ignored )
  ## - 20°C ( for this difference in temperature barely matters)
  # define Argon
  let ar = Argon.init()
  # and compute its density in CAST conditions
  let ρ_Ar = density(1050.mbar.to(Pascal), 293.K, ar.molarMass)
  # compute number density
  let n = numberDensity(ρ_Ar, ar.molarMass)
  result = absorptionLength(energy, n, ar.f2eval(energy)).to(cm)
