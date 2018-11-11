import jsffi, jsbind

# allows to convert `JsObject` to string
proc toString*(x: JsObject): cstring {.jsimportgWithName: "String".}
