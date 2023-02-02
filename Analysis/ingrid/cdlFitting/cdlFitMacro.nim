import macros, strutils, tables, math

import seqmath

type
  FitFuncKind* = enum
    ffConst, ffPol1, ffPol2, ffGauss, ffExpGauss

  FitFuncArgs* = object
    name*: string
    case kind*: FitFuncKind
    of ffConst:
      c*: float
    of ffPol1:
      cp*: float
    of ffPol2:
      cpp*: float
    of ffGauss:
      gN*: float
      gmu*: float
      gs*: float
    of ffExpGauss:
      ea*: float
      eb*: float
      eN*: float
      emu*: float
      es*: float


import macrocache
const FuncsTeX = CacheTable"TeXFunctions"
const fixed* = NaN

proc expGauss*(p: seq[float], x: float): float =
  # exponential * times (?!) from Christoph's expogaus.c
  if len(p) != 5:
    return Inf
  let p_val = 2.0 * (p[1] * pow(p[4], 2.0) - p[3])
  let q_val = 2.0 * pow(p[4], 2.0) * p[0] + pow(p[3], 2.0) - ln(p[2]) * 2.0 * pow(p[4], 2.0)

  let threshold = - p_val / 2.0 - sqrt( pow(p_val, 2.0) / 4.0 - q_val )

  if x < threshold:
    result = exp(p[0] + p[1] * x)
  else:
    result = p[2] * gauss(x, p[3], p[4])

proc genFitFuncImpl(idx: var int, xNode,
                    pFitNode: NimNode, paramsNode: seq[NimNode]): NimNode =
  ## the compile time procedure that creates the implementation lines for
  ## the <target><Filter>Funcs that we create, which calls the correct functions,
  ## e.g.
  ##   result += p_ar[i] * gauss(x, p[i + 1], p[i + 2])
  ##   inc i, 3
  let resultNode = ident"result"
  expectKind(pFitNode, nnkObjConstr)
  let fkind = parseEnum[FitFuncKind](pFitNode[2][1].strVal)
  case fKind
  of ffConst:
    let p0 = paramsNode[idx]
    result = quote do:
      `resultNode` += `p0`
    inc idx
  of ffPol1:
    let p0 = paramsNode[idx]
    result = quote do:
      `resultNode` += `p0` * `xNode`
    inc idx
  of ffPol2:
    let p0 = paramsNode[idx]
    result = quote do:
      `resultNode` += `p0` * `xNode` * `xNode`
    inc idx
  of ffGauss:
    let
      p0 = paramsNode[idx]
      p1 = paramsNode[idx + 1]
      p2 = paramsNode[idx + 2]
    result = quote do:
      `resultNode` += `p0` * gauss(`xNode`, `p1`, `p2`)
    inc idx, 3
  of ffExpGauss:
    let
      p0 = paramsNode[idx]
      p1 = paramsNode[idx + 1]
      p2 = paramsNode[idx + 2]
      p3 = paramsNode[idx + 3]
      p4 = paramsNode[idx + 4]
    result = quote do:
      `resultNode` += expGauss(@[`p0`, `p1`, `p2`, `p3`, `p4`], `xNode`)
    inc idx, 5

proc idOf(x: int): NimNode =
  let param = ident"p_ar"
  result = quote do:
    `param`[`x`]

proc drop(p: NimNode, frm: int): NimNode =
  result = copy(p)
  result.del(0, frm)
  #echo result.treeRepr

proc resolveCall(tab: OrderedTable[string, NimNode],
                 n: NimNode): NimNode =
  expectKind(n, nnkCall)
  const paramIdent = "p_ar"
  if n[1].strVal == paramIdent:
    # call is an nnkOpenSymChoice from the `p_ar[x]`
    result = n
  else:
    # call is a reference to a different parameter
    let key = n[0].strVal & "_" & n[1].strVal
    result = tab[key]

proc resolveNode(n: NimNode, tab: OrderedTable[string, NimNode]): NimNode =
  ## Recursively resolves all nodes in `n`
  result = copyNimTree(n)
  case n.kind
  of nnkCall:
    result = resolveCall(tab, n)
  else:
    for i in 0 ..< result.len:
      case result[i].kind
      of nnkIdent, nnkBracketExpr, nnkIntLit .. nnkFloatLit:
        # continue on all trivial types we like to keep as they are
        continue
      of nnkCall:
        # resolve calls (my contain e.g. `emu("Mn-Kalpha")` etc.
        result[i] = resolveCall(tab, result[i])
      of nnkInfix, nnkPar:
        # infix, recurse
        result[i] = resolveNode(result[i], tab)
      else:
        # else just break. Easier to add other functionality than allowing it and have it work
        # in a broken way
        error("Unsupported kind " & $result[i].kind & " with value: " & $result[i].repr)
  #case v.kind
  #of nnkCall:
  #  result.add resolveCall(tab, v)
  #else:
  #  if v.len == 0:
  #    result.add v
  #  else:
  #    # have to walk the tree to replace references
  #    var node = copyNimTree(v)
  #    echo node.treeRepr
  #    for i in 0 ..< node.len:
  #      echo "Node is ", node[i].kind, " : ", node[i].repr
  #      case node[i].kind
  #      of nnkCall:
  #        # replace the call with the resolved node
  #        node[i] = resolveCall(tab, node[i])
  #      else:
  #        # else leave node unchanged
  #        discard
  #    result.add node

proc resolveNodes(tab: OrderedTable[string, NimNode]): seq[NimNode] =
  ## resolve the potentially still ambiguous nodes (e.g. a reference to a
  ## parameter of a different line) and convert the `OrderedTable` into a
  ## serialized sequence of the nodes
  for k, v in tab:
    let ad = resolveNode(v, tab)
    result.add ad

proc incAndFill(tab: var OrderedTable[string, NimNode],
                idx: var int, keyBase, lineName: string) =
  let key = keyBase & "_" & lineName
  if key notin tab:
    tab[key] = idOf(idx)
    inc idx
  else:
    # if key in the table, remove it and re-add so that the position
    # is the correct one in our OrderedTable!
    let tmp = tab[key]
    tab.del(key)
    tab[key] = tmp

proc fillParamsTable(tab: var OrderedTable[string, NimNode],
                     p: NimNode, idx: var int) =
  ## creates a table of NimNodes, which stores the NimNode that will show up
  ## in the created fitting procedure under a unique string, which consists
  ## of `FitParameter_LineName`. This table can then later be used to
  ## reference arbitrary parameters of the fitting functions in a different
  ## or same line
  # get the correct part of the parameters
  let fkind = parseEnum[FitFuncKind](p[2][1].strVal)
  #var cTab = initOrderedTable[string, NimNode]()
  var toReplace: NimNode
  let lineName = p[1][1].strVal
  if p.len > 3:
    # if p.len has more than 3 children, one ore more parameters have
    # been fixed to some constant or relative value. Fill the `tab` with
    # the NimNodes corresponding to those values
    toReplace = p.drop(3)
    for x in toReplace:
      case x[1].kind
      of nnkIdent, nnkIntLit .. nnkFloatLit, nnkInfix, nnkCall:
        let id = x[0].strVal & "_" & lineName
        tab[id] = x[1]
      else: error("Unsupported kind " & $x[1].kind & " of val: " & $x[1].repr)

  # Now walk over the table with the parameters of the current line
  # and either enter the correct `p_ar[idx]` field as the value in
  # the table or put the value to be replaced at the correct position
  # in the table
  case fKind
  of ffConst:
    incAndFill(tab, idx, "c", lineName)
  of ffPol1:
    incAndFill(tab, idx, "cp", lineName)
  of ffPol2:
    incAndFill(tab, idx, "cpp", lineName)
    inc idx
  of ffGauss:
    incAndFill(tab, idx, "gN", lineName)
    incAndFill(tab, idx, "gmu", lineName)
    incAndFill(tab, idx, "gs", lineName)
  of ffExpGauss:
    incAndFill(tab, idx, "ea", lineName)
    incAndFill(tab, idx, "eb", lineName)
    incAndFill(tab, idx, "eN", lineName)
    incAndFill(tab, idx, "emu", lineName)
    incAndFill(tab, idx, "es", lineName)


proc serialize(n: NimNode): string

proc toElement(s: string): string =
  result = r"\ce{" & s & r"}"

proc toLine(s: string): string =
  ## Argument must be of the type `Kalpha`, `Lbeta` etc.
  result = $s[0]
  case s[1 .. ^1]
  of "alpha": result.add r"α"
  of "beta": result.add r"β"
  of "gamma": result.add r"γ"
  else:
    raise newException(ValueError, "Unknown line: " & $s)

proc toSubscript(s: string): string =
  result = r"_{" & s & r"}"

proc toEsc(s: string): string =
  case s
  of "esc": result = r"\text{esc}"
  else: doAssert false, "Invalid: " & $s

proc toSuperscript(s: string): string =
  result = r"^{" & s & r"}"

proc toLineName(s: string, attachTo = ""): string =
  if "-" notin s:
    result = r"_{\text{unknown}}"
  else:
    let parts = s.split("-")
    let element = toElement(parts[0])
    #if attachTo.len > 0:
    #  result = toSuperscript(parts[0]) #& attachTo # place element into a `^` then add
    #else:
    #  result = toSuperscript(parts[0])
    if parts.len == 2: # of type Cu-Kalpha
      let line = toLine(parts[1])
      if line == "esc":
        raise newException(ValueError, "Given line name is not valid. Missing the line (Kα, Lβ, ...) for argument: " & $s)
      result = toSuperscript(element) & toSubscript(line)
    elif parts.len == 3: # of type Cu-Kalpha-esc
      let line = toLine(parts[1])
      let esc = toEsc(parts[2])
      result = toSuperscript(element & ", " & esc) & toSubscript(line)

proc toFuncName(s: string): string =
  case s
  of "ffGauss": result = "G"
  of "ffExpGauss": result = "EG"

proc paramToTeX(s: string): string =
  case s
  of "a": result = "a"
  of "b": result = "b"
  of "mu": result = "μ"
  of "s": result = "σ"
  of "N": result = "N"

proc callToTex(n: NimNode): string =
  doAssert n[0].kind in {nnkSym, nnkIdent}
  let field = n[0].strVal
  let fn = field[0] # g or e
  doAssert fn in {'g', 'e'}
  var param = ""
  if fn == 'e': param = "e"
  param.add paramToTex(field[1 .. ^1])
  doAssert n[1].kind == nnkStrLit
  result = param & n[1].strVal.toLineName()

proc infixToTex(n: NimNode): string =
  doAssert n.kind == nnkInfix and n[0].kind in {nnkSym, nnkIdent} and n[0].strVal in ["/", "*"]
  case n[0].strVal
  of "/": result = r"\frac{" & n[1].serialize() & r"}{" & n[2].serialize() & r"}"
  of "*": result = n[1].serialize() & r"·" & n[2].serialize()
  else:
    doAssert false, "Unsupported infix kind: " & $n.repr

proc serialize(n: NimNode): string =
  case n.kind
  of nnkExprColonExpr:
    doAssert n[0].kind in {nnkSym, nnkIdent}
    case n[0].strVal
    of "name": result = n[1].strVal.toLineName() # name of the line
    of "kind": result = n[1].strVal.toFuncName() # which type of function it is
    else:
      result = n[1].serialize()
  of nnkInfix: result = n.infixToTex()
  of nnkCall: result = n.callToTex() # call must be a fake call, `gmu`, `gN`, `eN`, etc, reference to another line
  of nnkFloatLit: result = n.repr # repr is good here
  of nnkStrLit: doAssert false, "nnkStrLit should not appear on its own, only in `nnkCall`"
  of nnkSym, nnkIdent: doAssert false, "nnkIdent should not appear on its own"
  of nnkBracketExpr: result = n.repr # fine for this, should be `p_ar[i]`
  of nnkPar: result = r"(" & n[0].serialize() & r")"
  else: doAssert false, "Unsupported node kind: " & $n.kind & " in: " & $n.treerepr

proc toTex(p: NimNode): string =
  doAssert p.kind == nnkObjConstr and p[0].kind in {nnkIdent, nnkSym} and p[0].strVal == "FitFuncArgs", "No, was: " & $p.treerepr
  let lineName = p[1].serialize()
  let fnName = p[2].serialize()
  result = fnName & lineName
  if p.len > 3:
    result.add r"\left( "
  for i in 3 ..< p.len:
    result.add p[i].serialize()
    if i != p.len - 1:
      result.add ", "
  if p.len > 3:
    result.add r" \right)"

proc fnToTex(parts: NimNode): string =
  result = r"$"
  for i in 0 ..< parts.len:
    result.add parts[i].toTex()
    if i < parts.len - 1:
      result.add " + "
  result.add r"$"

proc buildFitFunc(name, parts: NimNode): NimNode =
  ## builds a CDL fit function based on the function described by
  ## the `seq[FitFuncArgs]` at compile time. Using the `FitFuncKind` of
  ## each part, it'll write the needed implementation lines for the
  ## call to the correct functions, e.g. `gauss`, `expGauss` etc.
  # define the variables needed in the implementation function and
  # for the parameters
  let
    idx = ident"i"
    paramsNode = ident"p_ar"
    xNode = ident"x"
  # define parameters and return type of the proc we create
  let
    retType = ident"float"
    retParNode = nnkIdentDefs.newTree(paramsNode,
                                      nnkBracketExpr.newTree(
                                        ident"seq",
                                        ident"float"),
                                      newEmptyNode())
    retXNode = nnkIdentDefs.newTree(xNode,
                                    ident"float",
                                    newEmptyNode())
  # create a node to hold the procedure body
  var procBody = newStmtList()

  var i = 0

  var paramsTab = initOrderedTable[string, NimNode]()
  for p in parts:
    # create the table which maps the parameter of each line to the
    # correct NimNode. May either be
    # - `p_ar[idx]`
    # - some literal int / float
    # - some calculation involving a reference to some other parameter
    paramsTab.fillParamsTable(p, i)
  # now that we have the whole tableresolve all still remaining references
  # to other parameters
  let params = resolveNodes(paramsTab)
  i = 0
  for p in parts:
    procBody.add genFitFuncImpl(i, xNode, p, params)

  # now define the result variable as a new proc
  result = newProc(name = nnkPostfix.newTree(ident"*", name),
                   params = [retType, retParNode, retXNode],
                   body = procBody,
                   procType = nnkFuncDef)


  FuncsTeX[name.strVal] = newLit(parts.fnToTex())

macro declareFitFunc*(name, stmts: untyped): untyped =
  ## DSL to declare the fit functions without having to go the
  ## const of `seq[FitFuncArgs]` + buildFitFunc macro route.
  ##
  ## .. code-block::
  ##   declareFitFunc(cEpic):
  ##     ffGauss: "C-Kalpha"
  ##     ffGauss: "O-Kalpha"
  ##   # will be expanded to
  ##   var cEpicF {.compileTime.} = @[FitFuncArgs(name: "C-Kalpha", kind: ffGauss),
  ##                                  FitFuncArgs(name: "O-Kalpha", kind: ffGauss)]
  ##   buildFitFunc(cEpicFunc, cEpicF)

  let
    thVarId = ident(name.strVal & "F_Mangle")
    funcId = ident(name.strVal & "Func")
  var ffSeq = nnkBracket.newTree()

  for s in stmts:
    expectKind(s, nnkCall)
    #echo s.treeRepr
    ## NOTE: while we could in principle use `FitFuncArgs` as a regular object at CT, we'd have
    ## to deal with the symbolic expressions for the fixed parameters somehow. That would make
    ## such a solution similarly complex to the AST one we currently use.
    let fkind = s[0]
    case s[1].len
    of 1:
      let ffName = s[1][0].strVal
      ffSeq.add quote do:
        FitFuncArgs(name: `ffName`, kind: `fKind`)
    else:
      var ffArg = nnkObjConstr.newTree(ident"FitFuncArgs")
      for x in s[1]:
        expectKind(x, nnkAsgn)
        if x[0].basename.ident == toNimIdent"name":
          ffArg.add nnkExprColonExpr.newTree(ident"name", x[1])
        else:
          ffArg.add nnkExprColonExpr.newTree(x[0], x[1])
      ffArg.insert(2, nnkExprColonExpr.newTree(ident"kind", fKind))
      ffSeq.add ffArg

  result = buildFitFunc(funcId, ffSeq)
  #echo result.repr

func filterNaN(s: openArray[float]): seq[float] =
  result = newSeqOfCap[float](s.len)
  for x in s:
    if classify(x) != fcNaN:
      result.add x

proc serialize*(parts: seq[FitFuncArgs]): seq[float] =
  for p in parts:
    case p.kind
    of ffConst:
      result.add filterNaN([p.c])
    of ffPol1:
      result.add filterNaN([p.cp])
    of ffPol2:
      result.add filterNaN([p.cpp])
    of ffGauss:
      result.add filterNaN([p.gN, p.gmu, p.gs])
    of ffExpGauss:
      result.add filterNaN([p.ea, p.eb, p.eN, p.emu, p.es])

proc toStr*(parts: seq[FitFuncArgs]): string =
  ## Generates a string

when isMainModule:
  import ingrid / calibration
  import ggplotnim, sequtils, seqmath
  # an example on how to define and use the above macro

  # start by declaring the fit function desired
  declareFitFunc(feSpec):
    ffExpGauss: "Mn-Kalpha-esc"
    ffExpGauss:
      name = "Mn-Kbeta-esc"
      eN = eN("Mn-Kalpha-esc") * p_ar[14]
      emu = emu("Mn-Kalpha-esc") * 3.5 / 2.9 # lock to relation to `Mn-Kalpha-esc` arg
      es = es("Mn-Kalpha-esc") # lock to `es` of `Mn-Kalpha`
    ffExpGauss: "Mn-Kalpha"
    ffExpGauss:
      name = "Mn-Kbeta"
      eN = eN("Mn-Kalpha") * p_ar[14]
      emu = emu("Mn-Kalpha") * 6.35 / 5.75 # lock to relation to `Mn-Kalpha` arg
      es = es("Mn-Kalpha") # lock to `es` of `Mn-Kalpha`

  #[ which will then generate the following function:
  func feSpecFunc(p_ar: seq[float]; x: float): float =
    result +=
        expGauss(@[p_ar[0], p_ar[1], p_ar[2], p_ar[3], p_ar[4]], x)
    result +=
        expGauss(@[p_ar[5], p_ar[6], p_ar[2] * p_ar[14], p_ar[3] * 3.5 / 2.9, p_ar[4]], x)
    result +=
        expGauss(@[p_ar[7], p_ar[8], p_ar[9], p_ar[10], p_ar[11]], x)
    result +=
        expGauss(@[p_ar[12], p_ar[13], p_ar[9] * p_ar[14], p_ar[10] * 6.35 / 5.75,
                   p_ar[11]], x)

  where we can see that all parameters of the first `expGauss` are still free and
  those of the second `expGauss` have either been fixed, are free or are tied to
  parameters of the first gauss.
  In order to fit or plot this function, we need some sensible parameters (or start
  parameters for a fit).
  For this we use the `FitFuncArgs` type. One element for each of the terms of the
  function.
  ]#
  let paramsFfa = @[
    FitFuncArgs(
      name: "Mn-Kalpha-esc",
      kind: ffExpGauss,
      ea: 1e-4,
      eb: 1e-5,
      eN: 500.0,
      emu: 95.0,
      es: 15.0),
    FitFuncArgs(
      name: "Mn-Kbeta-esc",
      kind: ffExpGauss,
      ea: 1e-4,
      eb: 1e-5,
      eN: fixed,
      emu: fixed,
      es: fixed), # additional parameters fixed, `fixed` is just an overload for `NaN`
    FitFuncArgs(
      name: "Mn-Kalpha",
      kind: ffExpGauss,
      ea: 1e-4,
      eb: 1e-5,
      eN: 2500.0,
      emu: 220.0,
      es: 15.0),
    FitFuncArgs(
      name: "Mn-Kbeta",
      kind: ffExpGauss,
      ea: 1e-4,
      eb: 1e-5,
      eN: fixed,
      emu: fixed,
      es: fixed) # additional parameters fixed
  ]
  # with this seq we can now generate a `seq[float]` by calling `serialize`:
  var params = paramsFfa.serialize
  #[ However, in the above declaration of our fitting function, we make use of
  an additional parameter 15, namely `p_ar[14]`. If we just serialize all `fixed`
  arguments will be removed (`params.len == 14` now). So we have to add an additional
  dummy parameter ]#
  # ratio of NAlpha/NBeta
  params.add 17.0 / 150.0
  # with this in hand we can now plot the function
  let x = linspace(0, 400, 1000)
  let y = x.mapIt(feSpecFunc(params, it))
  let df = toDf(x, y)
  ggplot(df, aes("x", "y")) +
    geom_line() +
    ggsave("feSpecFuncTest.pdf")
