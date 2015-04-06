################### u.: Phil's utility code
################### u.an: Phil's data analysis code
################### u.an.sr: Encapsulated file organization -- everything for StatRecord class

# To be sourced by u.an

################### StatRecord init ----

#### setClass: an.StatRecord
setClass (
  Class="an.StatRecord",
  representation(
    statName="character",
    dataNames="list",
    dataSource="character",
    desc="character",
    desc.set="character",
    ht="list"
  )
)

#### setClass: an.StatRecordSet
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

#### setClass: an.StatRecord.cor
setClass (
  Class="an.StatRecord.cor",
  representation(
    levene="list"
  ),
  contains="an.StatRecord"
)

#### setClass: an.StatRecord.reliability
setClass (
  Class="an.StatRecord.reliability",
  representation(
    stat.r.SBCorrected="numeric",
    cronbach.alpha="numeric"
  ),
  contains="an.StatRecord"
)

#### setClass: an.StatRecord.ttest
setClass (
  Class="an.StatRecord.ttest",
  representation(
    sd.pooled="numeric",
    levene="list"
  ),
  contains="an.StatRecord"
)
#### setClass: an.StatRecord.ttest.1sample
setClass (
  Class="an.StatRecord.ttest.1sample",
  representation(
    mu="numeric",
    sd="numeric"
  ),
  contains="an.StatRecord"
)

#### setClass: an.StatRecord.lme
setClass (
  Class="an.StatRecord.lme",
  representation(
  ),
  contains="an.StatRecord"
)

################### StatRecord generic methods ----
# Method: print()
setMethod ("print", signature="an.StatRecord",
           function (x) { # Note - the exact named param x is required
             print(getOutput(x), row.names=F)
           }
)
# Method: show()
setMethod ("show", signature="an.StatRecord",
           function (object) {
             return (print(object))
           }
)

# Method: getOutput() for an.StatRecord.cor
setMethod ("getOutput", signature="an.StatRecord.cor",
           function (x) {
             sr <- x
             ht <- an.sr.get.ht(sr)
             oe <- data.frame(vec.1=getSlot(sr, "dataNames")[[1]],
                              vec.2=getSlot(sr, "dataNames")[[2]],
                              stat.r=ht$estimate[[1]],
                              sig=ifelse(ht$p.value < PH$STAT.P.TRIPLESTAR.CUTOFF, "***", 
                                         ifelse(ht$p.value < PH$STAT.P.DOUBLESTAR.CUTOFF, "**",
                                                ifelse(ht$p.value < PH$STAT.P.SINGLESTAR.CUTOFF, "*",
                                                       ifelse(ht$p.value < PH$STAT.P.TREND.CUTOFF, "tnd", "_")))),
                              p=ht$p.value[[1]],
                              source=getSlot(sr, "dataSource"),
                              df=ht$parameter[[1]], # The numeric is to remove the name of the numeric, as that makes row names happen for some reason
                              desc=getSlot(sr, "desc"),
                              desc.set=getSlot(sr, "desc.set")
             )
             return (oe)
           }
)
# Method: getOutput() for an.StatRecord.reliability
# This is the same as for an.StatRecord.cor, except it adds in the SBCorrected column
setMethod ("getOutput", signature="an.StatRecord.reliability",
           function (x) {
             sr <- x
             ht <- an.sr.get.ht(sr)
             oe <- data.frame(vec.1=getSlot(sr, "dataNames")[[1]],
                              vec.2=getSlot(sr, "dataNames")[[2]],
                              p=ht$p.value[[1]],
                              r.halves=ht$estimate[[1]],
                              r.halves.SBCor=getSlot(sr, "stat.r.SBCorrected"),
                              cronb=getSlot(sr, "cronbach.alpha"),
                              source=getSlot(sr, "dataSource"),
                              df=ht$parameter[[1]], # The numeric is to remove the name of the numeric, as that makes row names happen for some reason
                              #                               sig=ifelse(ht$p.value < PH$STAT.P.TRIPLESTAR.CUTOFF, "***", 
                              #                                          ifelse(ht$p.value < PH$STAT.P.DOUBLESTAR.CUTOFF, "**",
                              #                                                 ifelse(ht$p.value < PH$STAT.P.SINGLESTAR.CUTOFF, "*",
                              #                                                        ifelse(ht$p.value < PH$STAT.P.TREND.CUTOFF, "tnd", "")))),
                              desc=getSlot(sr, "desc"),
                              desc.set=getSlot(sr, "desc.set")
             )
             return (oe)
           }
)
# Method: getOutput() for an.StatRecord.ttest
setMethod ("getOutput", signature="an.StatRecord.ttest",
           function (x) {
             sr <- x
             ht <- an.sr.get.ht(sr)
             p.levene <- getSlot(sr, "levene")[[1]][["Pr(>F)"]][1]
             oe <- data.frame( vec.1=do.call(paste, as.list(getSlot(sr, "dataNames")[[1]])),
                               vec.2=do.call(paste, as.list(getSlot(sr, "dataNames")[[2]])),
                               sig=pa(ifelse(ht$p.value < PH$STAT.P.TRIPLESTAR.CUTOFF, "***", 
                                             ifelse(ht$p.value < PH$STAT.P.DOUBLESTAR.CUTOFF, "**",
                                                    ifelse(ht$p.value < PH$STAT.P.SINGLESTAR.CUTOFF, "*",
                                                           ifelse(ht$p.value < PH$STAT.P.TREND.CUTOFF, "tnd", "")))),
                                      ifelse(p.levene < PH$STAT.LEV.CUTOFF, "LEV", "")),
                               p=ht$p.value[[1]],
                               source=getSlot(sr, "dataSource"),
                               M.diff=ht$estimate[["mean of x"]] - ht$estimate[["mean of y"]],
                               M.1=ht$estimate[["mean of x"]],
                               M.2=ht$estimate[["mean of y"]],
                               SDpool=getSlot(sr, "sd.pooled"),
                               d=(ht$estimate[["mean of x"]] - ht$estimate[["mean of y"]]) / getSlot(sr, "sd.pooled"),
                               CI.name.L=ht$conf.int[1],
                               CI.name.U=ht$conf.int[2],
                               p.levene=p.levene,
                               df=ht$parameter[[1]],
                               stat.t=ht$statistic[[1]],
                               desc=getSlot(sr, "desc"),
                               desc.set=getSlot(sr, "desc.set")
             )
             CI.name.L.new <- pa("CI", attributes(ht$conf.int)[["conf.level"]], "L")
             CI.name.L.new <- strs.replace(CI.name.L.new, "0.", ".")
             names(oe)[names(oe)=="CI.name.L"] <- CI.name.L.new
             names(oe)[names(oe)=="CI.name.U"] <- strs.replace(CI.name.L.new, "L", "U")
             return (oe)
           }
)

# Method: getOutput() for an.StatRecord.ttest
setMethod ("getOutput", signature="an.StatRecord.ttest.1sample",
           function (x) {
             sr <- x
             ht <- an.sr.get.ht(sr)
             oe <- data.frame( vec=getSlot(sr, "dataSource"),
                               sig=pa(ifelse(ht$p.value < PH$STAT.P.TRIPLESTAR.CUTOFF, "***", 
                                             ifelse(ht$p.value < PH$STAT.P.DOUBLESTAR.CUTOFF, "**",
                                                    ifelse(ht$p.value < PH$STAT.P.SINGLESTAR.CUTOFF, "*",
                                                           ifelse(ht$p.value < PH$STAT.P.TREND.CUTOFF, "tnd", ""))))),
                               p=ht$p.value[[1]],
                               mean=ht$estimate[["mean of x"]],
                               mu=getSlot(sr, "mu"),
                               d=(ht$estimate[["mean of x"]] - getSlot(sr, "mu")) / getSlot(sr, "sd"),
                               SD=getSlot(sr, "sd"),
                               CI.name.L=ht$conf.int[1],
                               CI.name.U=ht$conf.int[2],
                               source=do.call(paste, as.list(getSlot(sr, "dataNames")[[1]])),
                               df=ht$parameter[[1]],
                               stat.t=ht$statistic[[1]],
                               desc=getSlot(sr, "desc"),
                               desc.set=getSlot(sr, "desc.set")
             )
             CI.name.L.new <- pa("CI", attributes(ht$conf.int)[["conf.level"]], "L")
             CI.name.L.new <- strs.replace(CI.name.L.new, "0.", ".")
             names(oe)[names(oe)=="CI.name.L"] <- CI.name.L.new
             names(oe)[names(oe)=="CI.name.U"] <- strs.replace(CI.name.L.new, "L", "U")
             return (oe)
           }
)

################### StatRecord get methods ----
# Method: an.sr.get.ht()
# Use this to get htest from ht slot instead of list(htest)
setGeneric ("an.sr.get.ht",
            def=function (sr) { standardGeneric("an.sr.get.ht") })
setMethod ("an.sr.get.ht", signature="an.StatRecord",
           function (sr) {
             return (getSlot(sr, "ht")[[1]])
           }
)

################### StatRecord creation methods ----

# Creates correlation StatRecord
anpriv.sr.create.cor <- function (dataSource="", dataNames="", desc="", desc.set="", ht=list(), levene=list()) {
  return (new(Class="an.StatRecord.cor",
              statName="cor",
              dataNames=dataNames,
              dataSource=dataSource,
              desc=desc,
              desc.set=desc.set,
              ht=ht,
              levene=levene))
}

# Creates reliability StatRecord
anpriv.sr.create.reliability <- function (dataSource="", dataNames="", desc="", desc.set="", stat.r.SBCorrected=NA, cronbach.alpha, ht=list()) {
  return (new(Class="an.StatRecord.reliability",
              statName="cor",
              dataNames=dataNames,
              dataSource=dataSource,
              desc=desc,
              desc.set=desc.set,
              ht=ht,
              stat.r.SBCorrected=stat.r.SBCorrected,
              cronbach.alpha=cronbach.alpha))
}

# Creates ttest, independent samples, StatRecord
anpriv.sr.create.ttest <- function (dataSource="", dataNames="", desc="", desc.set="", ht=list(), sd.pooled=as.numeric(NA), levene=list()) {
  return (new(Class="an.StatRecord.ttest",
              statName="ttest",
              dataNames=dataNames,
              dataSource=dataSource,
              desc=desc,
              desc.set=desc.set,
              ht=ht,
              sd.pooled=sd.pooled,
              levene=levene))
}

# Creates ttest, 1sample, StatRecord
anpriv.sr.create.ttest.1sample <- function (dataSource="", dataNames="", desc="", desc.set="", ht=list(), mu=0, sd=NA) {
  return (new(Class="an.StatRecord.ttest.1sample",
              statName="ttest.1s",
              dataNames=dataNames,
              dataSource=dataSource,
              desc=desc,
              desc.set=desc.set,
              ht=ht,
              mu=mu,
              sd=sd))
}

# Creates StatRecordSet, to be added on to
an.srSet.create <- function (desc.set="") {
  return (new(Class="an.StatRecordSet",
              statRecords=list(),
              desc.set=desc.set))
}

# end ----