################### u.: Phil's utility code
################### u.an: Phil's data analysis code
################### u.an.srSet: Encapsulated file organization -- everything for StatRecord class

################### StatRecordSet init ----
setClass (
  Class="an.StatRecordSet",
  representation(
    statRecords="list",
    desc.set="character"
  ),
  validity=function (object) {
    if (length(object@statRecords) > 0) {
      for (statRecord in object@statRecords) {
        if (class(statRecord) != class(object@statRecords[[1]]))
          stop("an.StatRecordSet validity check failed: statRecords are of multiple types")
        if (!inherits(statRecord))
          stop("an.StatRecordSet validity check failed: statRecords are not of class (inherited) an.StatRecord")
      }
    }
  }
)

################### StatRecord generic methods ----
# Prints StatRecordSet
setMethod ("print", signature="an.StatRecordSet",
           function (x) { # Note - the exact named param x is required
             print(getOutput(x), row.names=F)
           }
)
setMethod ("show", signature="an.StatRecordSet",
           function (object) {
             return (print(object))
           }
)

setMethod ("getOutput", signature="an.StatRecordSet",
           function (x) {
             srSet <- x
             oe <- data.frame()
             for (sr in srSet@statRecords) {
               oe <- rbind(oe, getOutput(sr))
             }
             return (oe)
           }
)

################### StatRecord add. methods ----
# Adds a StatRecord to a StatRecordSet
setGeneric ("an.srSet.add",
            def=function (srSet, sr) { standardGeneric("an.srSet.add") })
setMethod ("an.srSet.add", signature="an.StatRecordSet",
           function (srSet, sr) {
             srSet@statRecords <- c(srSet@statRecords, list(sr))
             return (srSet)
           }
)

################### StatRecord creation methods ----

# Creates StatRecordSet, to be added on to
an.srSet.create <- function (desc.set="") {
  return (new(Class="an.StatRecordSet",
              statRecords=list(),
              desc.set=desc.set))
}

