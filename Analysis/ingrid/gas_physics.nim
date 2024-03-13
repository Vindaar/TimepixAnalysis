import pkg / [xrayAttenuation, unchained, datamancer]

proc initCASTGasMixture*(): GasMixture =
  ## Returns the absorption length for the given energy in keV for CAST
  ## gas conditions:
  ## - Argon / Isobutane 97.7 / 2.3 %
  ## - 20°C ( for this difference in temperature barely matters)
  # define Argon
  let arC = compound((Ar, 1)) # need Argon gas as a Compound
  let isobutane = compound((C, 4), (H, 10))
  # define the gas mixture
  result = initGasMixture(293.K, 1050.mbar, [(arC, 0.977), (isobutane, 0.023)])

proc absorptionLengthCAST*(energy: keV): cm =
  ## Returns the absorption length for the given energy in keV for CAST
  ## gas conditions:
  ## - Argon / Isobutane 97.7 / 2.3 %
  ## - 20°C ( for this difference in temperature barely matters)
  # define Argon
  let arC = compound((Ar, 1)) # need Argon gas as a Compound
  let isobutane = compound((C, 4), (H, 10))
  # define the gas mixture
  let gm = initGasMixture(293.K, 1050.mbar, [(arC, 0.977), (isobutane, 0.023)])
  # and compute λ
  result = absorptionLength(gm, energy).to(cm)

proc absorptionLengthCAST*(gm: GasMixture, energy: keV): cm =
  ## Returns the absorption length for the given energy in keV for CAST
  ## gas conditions:
  ## - Argon / Isobutane 97.7 / 2.3 %
  ## - 20°C ( for this difference in temperature barely matters)
  result = absorptionLength(gm, energy).to(cm)
