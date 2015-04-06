################### u.: Phil's utility code
################### u.an: Phil's data analysis code
################### u.an.meth: User-accessible methods
################### u.an.meth Master: User-accessible methods code files

# To be sourced by u.an.R Master



# Method: getSlot()
# Generic accessor method, for any object
# The point of this is that if you want to change what happens upon accessing in the future, you can.
setGeneric ("getSlot",
            def=function (object, slotName) { standardGeneric("getSlot") })
setMethod ("getSlot", signature="ANY",
           function (object, slotName) {
             return (slot(object, slotName))
           })

# Method: getOutput()
# Generic method to get the output intended for any object... basically should give you what print() would give you, but without printing it.
setGeneric ("getOutput",
            def=function (x) { standardGeneric("getOutput") })
setMethod ("getOutput", signature="ANY",
           function (x) {
             return (x)
           })

setGeneric ("getOutput",
            def=function (x) { standardGeneric("getOutput") })
setMethod ("getOutput", signature="list",
           function (x) {
             return (lapply(x, pa, "\n"))
           })

setMethod ("getOutput", signature="data.frame",
           function (x) {
             return (x)
           })
