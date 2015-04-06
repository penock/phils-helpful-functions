################### PH: Phil's Helpful (functions)
# Instructions for use:
#   - Set the constant PATH.PH below to the correct path containing this file
#   - Source this file to get access to the functions
#
# Navigating this source file (with RStudio):
#   - My strategy has been to pack everything into this one file
#   - I did so because RStudio makes it easy to manage and navigate in the one file (or multiple, I guess, but I chose to do one file)
#   - Use keyboard shortcut for RStudio's Code menu -> "Jump to" to navigate to them (Alt-Shift-J), seeing both sections and function names
#   - Also, leave this codefile open when working on your projects. That way, you can press F2 while cursor is on a function name in your code,
#     and it'll take you here. Code -> Go to File/Function is useful too (Ctrl-.), as is typing ph. and Tab to see functions listed
PATH.PH <- "C:/Dropbox/PH/ph.src"



# ------------------------------------ < User options > ------------------------------------

PH <- list() # About the namespace:
             # The PH list holds all constants for PH, so that your variable namespace doesn't get crowded with these constants.

PH$REMOVE.PH.PREFIX <- TRUE # All ph functions start with "ph.", but when you set TRUE here, you can then access them without typing ph. in front. See function ph.removePHPrefix().

PH$SPEC.VERSION <- 0.2 # I try to make my functions backward-compatible. But sometimes, I need to change them.
                       # For example, in version 0.1, say my t-tests (stat.t) outputted 7 columns. If, in 0.2, they output 9 columns, well,
                       # if this would mess up your code, you can put in spec version 0.1, and where possible, it'll use the old specification.

PH$DATAFRAME.DISPLAY.WIDE.ENABLE <- TRUE  # To understand this, see below, the print.data.frame() function... that replaces the default, so it prints wider,
PH$DATAFRAME.DISPLAY.WIDE.WIDTH  <- 600   # and you just scroll right to view the data frame, much better than how they wrap by default in my opinion
PH$DATAFRAME.DISPLAY.ROW.NAMES   <- FALSE # This overrides the default R habit of putting row names (numbers) down the left side by making print.data.frame(row.names=FALSE) as default instead of TRUE.

PH$STAT.P.TREND.CUTOFF <- .1        # In stat output, under "sig" column, mark "tnd" for any p values under .1
PH$STAT.P.SINGLESTAR.CUTOFF <- .05  # ... * for any p < .05
PH$STAT.P.DOUBLESTAR.CUTOFF <- .01  # ... ** for any p < .01
PH$STAT.P.TRIPLESTAR.CUTOFF <- .001 # ... *** for any p < .001
PH$STAT.R.SINGLESTAR.CUTOFF <- .3   # In stat.cor output, mark 
PH$STAT.R.DOUBLESTAR.CUTOFF <- .6   # ... 
PH$STAT.R.TRIPLESTAR.CUTOFF <- .9   # ... 
PH$STAT.LEVENE.CUTOFF <- .05        # In stat.t output, under "sig" column, mark "LEV" if levene's test showed differing variances (though one generally ignores levene if sample large enough)

PH$FORMAT.DECIMALS <- 4 # 3 or 4 digits are both reasonable choices for most purposes

options(digits=PH$FORMAT.DECIMALS)
options(error = recover)

if (PH$DATAFRAME.DISPLAY.WIDE.ENABLE | PH$DATAFRAME.DISPLAY.ROW.NAMES) {
  print.data.frame <- function(..., row.names=PH$DATAFRAME.DISPLAY.ROW.NAMES, wide=T) {
    p <- get("print.data.frame", envir=baseenv())
    if (wide==T & options("width") < PH$DATAFRAME.DISPLAY.WIDE.WIDTH)
      options(width=PH$DATAFRAME.DISPLAY.WIDE.WIDTH)
    p(..., row.names=row.names)
  }
}

# ------------------------------------ < Autoexec and Functions to manage PH functions > -----------------------------

# EnviHack to get the function names etc. out of the Environment pane of RStudio (because causes RStudio it to lag),
# but keep them attached so they work just the same.
# I should be very wary of this, though, since attach() can make crazy things happen.
# I just couldn't tolerate the lag, and switching tabs away from Environment pane doesn't work because it regains focus upon browser() and recover().
PH$LS.BEFORE <- ls()


# Makes a copy of all functions with the passed prefix to without the prefix
# So, instead of typing ph.op(df), you can just type op(df)
# The exceptions are the following, which never have ph. in front of them (I use them too frequently to always type ph. when coding in this file):
# ca, cal, pa, pal, %_%, opm and functions beginning with opm
ph.removePHPrefix <- function (prefixToRemove = "ph.") {
  allObjects <- ls(envir=globalenv())
  toCopy <- allObjects[ph.strs.detect.prefix.vec(allObjects, 'ph.')]
  sapply(toCopy, FUN=function(funName) {
    newName <- ph.strs.replace(funName, 'ph.', '')
    assign(newName, get(funName), envir=globalenv())
  })
  invisible(NULL)
}

# ------------------------------------ < Load required libraries > ------------------------------------

# Note to user -- you can disable some of these loads if you want, if you're not using PH functions that need them
library(reshape)
library(stringr)
library(prettyR)
library(plyr)
library(WRS)
library(car)


if (!exists("PH$DISABLE.SOME.OUTPUT")) PH$DISABLE.SOME.OUTPUT <- F
assign("opm.lastopmsTime", Sys.time(), envir=baseenv())


# Returns vector of names of every function present in the file, based on whatever precedes the <- then "function" on any lines
# (Includes any named functions declared inside other functions.)
ph.sourceFile.listFunctions <- function(sourceFile, functionIdentifier = ' <- '%_%'function') {
  maxFunctionNameLength <- 100
  sourceText <- ph.read.entireFile(sourceFile)
  searchRes <- gregexpr(functionIdentifier, sourceText, fixed=T)
  hitLocs <- searchRes[[1]]
  
  # On each location, extract the function name as whatever comes before the function identifier
  functionNames <- sapply(hitLocs, function (hitLoc, sourceText) {
    preHitText <- substr(sourceText, hitLoc - maxFunctionNameLength, hitLoc-1)
    newlineSearchRes <- gregexpr('\n', preHitText, fixed=T)
    newlineLocs <- newlineSearchRes[[1]]
    startOfName <- newlineLocs[length(newlineLocs)] + 1
    substr(preHitText, startOfName, nchar(preHitText)) # returns name of function
  }, sourceText=sourceText)
  functionNames
}

# ------------------------------------ < Stat functions > ------------------------------------

source(file.path(PATH.PH, "/u.an.meth.R"))
source(file.path(PATH.PH, "/u.an.sr.R"))
source(file.path(PATH.PH, "/u.an.srSet.R"))


# Performs all possible numeric ttests of each vector comprising the 3rd dimension of each 2D cell on a passed list of oes.
# If you consider the oes to be all stacked up on the floor, each vector is like a vertical tube coming up in the ground.
# The oes must all match in shape and type of content.
# The first column must contain what you want to serve as the row names.
# Returns: a StatRecordSet containing t-test StatRecords
ph.oea.3D.ttests <- function (oeList, desc.set="3D ttests") {
  oe <- ph.oea.3D.to2D(oeList)
  oe <- ph.oe.remove.cols(oe, remove.ifAllSame=T, keepOnly.ifType="numeric")
  return (ph.oea.t.1sample.set(data=oe, colNames=names(oe), desc.set=desc.set))
}

# Takes a list of same-format oes and turns it into a data frame, where each column comes from a cell's 3rd dimension
# If you consider the oes to be all stacked up on the floor, each column in the output oe now represents a vertical tube coming up from the ground in the original stack.
# So, whereas the stack was original a stack of oes of dimensions X and Y, the output oe is ONE oe of dimension Z as rows and dimension X and Y labels unfolded into names(output oe)
# Intended for use by ph.an.oe.3D.ttests()
ph.oea.3D.to2D <- function (oeList) {
  oeList <- ph.list.list.promoteSingleElementLists(oeList)
  if (length(oeList) < 2)
    stop("ph.an.oe.3D.ttests(): You passed in less than 2 oes.")
  for (i in 1:length(oeList)) {
    ioe <- oeList[[i]]
    col1 <- ioe[1]
    restoe <- ioe[2:ncol(ioe)]
    restoe <- ph.oe.remove.cols(restoe, keepOnly.ifType="numeric")
    oeList[[i]] <- cbind(col1, restoe)
  }
  oe1 <- oeList[[1]]
  for (ioe in oeList) {
    if (nrow(ioe) != nrow(oe1) | ncol(ioe) != ncol(oe1))
      stop("ph.an.oe.3D.ttests(): oes are not of equal dimensions.")
    if (!identical(names(ioe), names(oe1)))
      warning("ph.an.oe.3D.ttests() WARNING: oes do not have identical column names.")
    if (!identical(ioe[[1]], oe1[[1]]))
      stop("ph.an.oe.3D.ttests(): first columns of oes are not identical. They should be what you want as row names, so should be identical.")
  }
  vecList <- list()
  for (ioe in oeList) {
    for (iColName in names(ioe)[-1]) {
      for (iRowNum in 1:nrow(ioe)) {
        cellName <- pa(iColName, "__", ioe[[iRowNum, 1]])
        if (!(cellName %in% names(vecList)))
          vecList[cellName] <- ioe[[iRowNum, iColName]]
        else
          vecList[[cellName]] <- c(vecList[[cellName]], ioe[[iRowNum, iColName]])
      }
    }
  }
  return (as.data.frame(vecList))
}


# My correlation SET
ph.oea.cor.set <- function (data, colNames.left, colNames.top=colNames.left, desc.set="cor set", dataSource=as.character(match.call()[[2]]), noDoubles=T, ...) {
  if (length(ph.listAbsences(names(data), colNames.left)) > 0)
    stop(pa("ph.an.cor.set(): these colNames you asked to test don't actually exist: ", pa(ph.listAbsences(names(data), colNames.left), collapse=", ")))
  if (length(ph.listAbsences(names(data), colNames.top)) > 0)
    stop(pa("ph.an.cor.set(): these colNames you asked to test don't actually exist: ", pa(ph.listAbsences(names(data), colNames.top), collapse=", ")))
  srSet <- an.srSet.create(desc.set=desc.set)
  donePairs <- list()
  for (iColName.left in colNames.left) {
    for (iColName.top in colNames.top) {
      if (iColName.left == iColName.top)
        next
      if (noDoubles==T &
            any(sapply(donePairs, function (pair) identical(pair, c(iColName.top, iColName.left)))))
        next
      if (ph.vec.is.goodForTests(data[[iColName.left]])==F | 
            ph.vec.is.goodForTests(data[[iColName.top]])==F | 
            identical(data[[iColName.left]], data[[iColName.top]]))
        next
      donePairs[[length(donePairs)+1]] <- c(iColName.left, iColName.top)
      sr <- ph.stat.cor(data[[iColName.left]],
                   data[[iColName.top]],
                   srcNames=list(iColName.left, iColName.top),
                   dataSource=dataSource,
                   desc.set=desc.set,
                   ...)
      srSet <- an.srSet.add(srSet, sr)
    }
  }
  return (getOutput(srSet))
}


# My cohen's d, for a set of one group's pre-posts
# NEEDS FIXING: You must do na.omit across pre and post vars
# oea.coh.trt.set.btgrps<-function (data, varList, preSuffix, postSuffix, cond1, cond2) {
#   for (iVar in varList) {
#     preTrtvec.1 <- data[cond1, pa(iVar, preSuffix)]
#     postTrtvec.1 <- data[cond1, pa(iVar, postSuffix)]
#     cohens1 <- (mean(postTrtvec.1, na.rm=T) - mean(preTrtvec.1, na.rm=T)) / sd(preTrtvec.1, na.rm=T)
#     preTrtvec.2 <- data[cond2, pa(iVar, preSuffix)]
#     postTrtvec.2 <- data[cond2, pa(iVar, postSuffix)]
#     cohens2 <- (mean(postTrtvec.2, na.rm=T) - mean(preTrtvec.2, na.rm=T)) / sd(preTrtvec.2, na.rm=T)
#     ca("\nd diff=", format(cohens1 - cohens2, digits=3), " from group1 (n=", length(preTrtvec.1), ", d=", format(cohens1, digits=3), ") - group2 (n=3", length(preTrtvec.2), ", d=", format(cohens2, digits=3), ") : ", iVar, " ", postSuffix, "-", preSuffix)
#   }
# }
# an.coh.trt.s.btgrps<-function (...) oea.coh.trt.set.btgrps(...) # Legacy

# My cohen's d, for a set of one group's pre-posts
ph.oea.coh.trt.set <- function (data, stems, preSuffix, postSuffix) {
  ph.outoe.make(c("Stem", "d", "M."%_%preSuffix, "M."%_%postSuffix, "M.Diff", "SD."%_%preSuffix))
  
  for (stem in stems) {
    oe.forOmit <- na.omit(data[c(pa(stem, preSuffix), pa(stem, postSuffix))])
    preTrtVec <- oe.forOmit[[pa(stem, preSuffix)]]
    postTrtVec <- oe.forOmit[[pa(stem, postSuffix)]]
    cohens <- (mean(postTrtVec) - mean(preTrtVec)) / sd(preTrtVec)
    ph.outoe.push(row=data.frame(
      stem,
      cohens,
      mean(preTrtVec),
      mean(postTrtVec),
      mean(postTrtVec) - mean(preTrtVec),
      sd(preTrtVec)))
  }
  return (ph.outoe.get())
}


# My t-test, independent samples, SET
ph.oea.t.set <- function (data, colNames, grpFac=NULL, condLogiVec.1=NULL, condLogiVec.2=NULL, srcNames=NULL, desc.set="t set", smartSkip=T, na.rm=T, ...) {
  if (!is.null(grpFac)) {
    if (!is.factor(data[[grpFac]]) | length(levels(data[[grpFac]])) != 2)
      data[[grpFac]] <- factor(data[[grpFac]])
    if (length(levels(data[[grpFac]])) != 2)
      stop("ph.oea.t.set(): You need exactly 2 levels in your grpFac, but you actually have: "%_%length(levels(data[[grpFac]])))
    condLogiVec.1 <- data[[grpFac]] == levels(data[[grpFac]])[1]
    condLogiVec.2 <- data[[grpFac]] == levels(data[[grpFac]])[2]
    if (is.null(srcNames)) {
      srcNames <- list(levels(data[[grpFac]])[1], levels(data[[grpFac]])[2])
    }
  }
  if (length(ph.listAbsences(names(data), colNames)) > 0)
    stop(pa("ph.an.t.set(): these colNames you asked to test don't actually exist: ", pa(ph.listAbsences(names(data), colNames), collapse=", ")))
  if (is.na(srcNames[1]))
    srcNames <- list(ph.str.stripSpaces(deparse(match.call()[[4]])), ph.str.stripSpaces(deparse(match.call()[[5]])))
  srSet <- an.srSet.create(desc.set=desc.set)
  for (iColName in colNames) {
    if (ph.vec.is.goodForTests(data[condLogiVec.1, iColName])==F | 
          ph.vec.is.goodForTests(data[condLogiVec.2, iColName])==F |
          identical(data[condLogiVec.1, iColName], data[condLogiVec.2, iColName]))
      next
    sr <- ph.stat.t(data[condLogiVec.1, iColName],
               data[condLogiVec.2, iColName],
               srcNames=srcNames,
               desc=iColName,
               desc.set=desc.set,
               na.rm=na.rm,
               ...
    )
    srSet <- an.srSet.add(srSet, sr)
  }
  return (getOutput(srSet))
}

# My one-sample t-test SET
ph.oea.t.1sample.set <- function (data, colNames, desc.set="1s t set", ...) {
  srSet <- an.srSet.create(desc.set=desc.set)
  if (length(ph.listAbsences(names(data), colNames)) > 0)
    stop(pa("ph.an.t.1sample.set(): these colNames you asked to test don't actually exist: ", pa(ph.listAbsences(names(data), colNames), collapse=", ")))
  for (iColName in colNames) {
    if (ph.vec.is.goodForTests(data[[iColName]])==F)
      next
    sr <- ph.stat.t.1sample(vec=data[[iColName]],
                       dataNames=list(ph.str.stripSpaces(deparse(match.call()[[2]]))),
                       dataSource=iColName,
                       desc.set=desc.set,
                       ...)
    srSet <- an.srSet.add(srSet, sr)
  }
  return (getOutput(srSet))
}


# I set default to "the original levene's test" (by doing center=mean) because that's what SPSS uses
ph.an.leveneTest.2vecs <- function (vec.1, vec.2, center=mean) {
  oe <- data.frame(vecs=vec.1, group=1)
  oe <- rbind(oe, data.frame(vecs=vec.2, group=2))
  return (leveneTest(y=oe$vecs, group=as.factor(oe$group), center=mean))
}



# My cohen's d, for one group's pre-posts at a time
#   designed for treatment stats... it's standardized mean gain, but uses pre-treatment SD not pooled SD, as per what Rich says
ph.an.coh.trt <- function (preTrtVec, postTrtVec) {
  cohens <- (mean(postTrtVec, na.rm=T) - mean(preTrtVec, na.rm=T)) / sd(preTrtVec, na.rm=T)
  ca("d=", format(cohens, digits=3))
}



# My correlation, simplified with one-line df output
ph.stat.cor <- function (src.1, src.2, ..., desc='', note='', srcNames=c(ph.str.stripSpaces(deparse(match.call()[[2]])), ph.str.stripSpaces(deparse(match.call()[[3]])))) {
  res <- cor.test(src.1, src.2, ...)
  rowoe <- data.frame(src.1=srcNames[[1]],
                      src.2=srcNames[[2]],
                      stat.r=res$estimate[[1]],
                      sig=ifelse(res$p.value < PH$STAT.P.TRIPLESTAR.CUTOFF, '***', 
                                 ifelse(res$p.value < PH$STAT.P.DOUBLESTAR.CUTOFF, '**',
                                        ifelse(res$p.value < PH$STAT.P.SINGLESTAR.CUTOFF, '*',
                                               ifelse(res$p.value < PH$STAT.P.TREND.CUTOFF, 'tnd', '_')))),
                      p=res$p.value,
                      df=res$parameter[[1]],
                      desc=desc,
                      note=note)
  rowoe
}

# LEGACY: For the SHR reliability analysis code to work... an.cor() gets called by an.reliability.SHR() within u.an.func.R
ph.an.cor <- function (vec.1, vec.2, ..., desc="", desc.set="", dataSource="", dataNames=list(str.stripSpaces(deparse(match.call()[[2]])), str.stripSpaces(deparse(match.call()[[3]])))) {
  # Used to have a param na.rm and this... but I think we're better off letting cor.test do its default na.rm behavior, which appears to be pairwise deletion
  #   if (na.rm==T) {
  #     vec.1 <- vec.1[!is.na(vec.1)]
  #     vec.2 <- vec.2[!is.na(vec.2)]
  #   }
  res <- cor.test(vec.1, vec.2, ...)
  sr <- anpriv.sr.create.cor(dataNames=dataNames,
                             dataSource=dataSource,
                             desc=desc,
                             desc.set=desc.set,
                             ht=list(res),
                             levene=list(an.leveneTest.2vecs(vec.1, vec.2))
  )
  return (sr)
}

# Spearman-Brown Prophecy Formula
# Source: http://en.wikipedia.org/wiki/Spearman%E2%80%93Brown_prediction_formula
ph.an.formula.spearmanbrown.prophecy <- function (r.current, augmentation.ratio=2) {
  r <- r.current
  N <- augmentation.ratio
  return (
    (N*r) / (1 + (N - 1)*r)
  )
}

# Cronbach's Alpha formula
# Source: http://en.wikipedia.org/wiki/Cronbach's_alpha (where it says "Alternatively, the Cronbach's  can also be defined as")
# Note - I confirmed that this result is the same as the alpha() function in CRAN's psych package, the raw_alpha result, same I checked to 9 decimal places.
ph.an.formula.cronbach.alpha <- function (oe) {
  covmat <- cov(oe)
  # Iterate through half (a diagonal half) of the covariance matrix,
  # extracting the diagonals as variances and the non-diagonals as covariances
  item.variances <- vector()
  item.covariances <- vector()
  for (iRow in 1:nrow(covmat)) {
    for (iCol in (iRow:ncol(covmat))) {
      if (iRow == iCol)
        item.variances <- c(item.variances, covmat[iRow, iCol])
      else
        item.covariances <- c(item.covariances, covmat[iRow, iCol])
    }
  }
  
  k <- ncol(oe)
  c <- mean(item.covariances)
  v <- mean(item.variances)
  alpha <- (k * c) / (v + (k-1) * c)
  return (alpha)
}


# Missing data handling... fills in LOCF in the oe and returns the modified oe
ph.an.mis.fill.LOCF <- function (oe, exts, stems="[ALL THAT HAVE FIRST EXT]", stepName="1", silent=F) {
  # To get "[ALL]" stems, grab any colname that ends with exts[1]
  if (stems[1] == "[ALL THAT HAVE FIRST EXT]") {
    stems.full <- names(oe)[ph.strs.detect.suffix.vec(names(oe), exts[1])]
    stems <- strs.replace.vec(stems.full, exts[1], "")
    if (silent==F)
      opm("For LOCF filling, found stems: "%_%pa(stems, collapse=", "))
  }
  
  if (length(exts) < 2)
    stop("ph.an.mis.fill.LOCF(): must have at least 2 exts")
  
  oe[[pa("LOCF.count.step.", stepName)]] <- 0 # This col represents the number of fill-ins this row has had
  for (stem in stems) {
    # Make sure the stem's columns don't skip any timepoints... but it's ok if exts are like w0 through f2, and only w0 through w4 are valid for a stem
    valid <- ph.oe.ext.validExts(oe, stem, exts)
    for (i in seq_along(valid)) {
      if (valid[i] != exts[i])
        stop("ph.an.mis.fill.LOCF(): looks like ext '"%_%exts[i]%_%"' doesn't exist for stem '"%_%stem%_%"'")
    }
    # Increment up the exts, for this stem, filling them in if missing. (Leaves them NA if points from current ext back to first ext are all NA.)
    for (i in 2:length(exts)) {
      curColName <- stem%_%exts[i]
      curColVec.old <- oe[[curColName]]
      oe[is.na(curColVec.old), curColName] <- oe[is.na(curColVec.old), stem%_%exts[i-1]]
      curColVec.filled <- oe[[curColName]]
      oe[[pa("LOCF.count.step.", stepName)]][is.na(curColVec.old) & !is.na(curColVec.filled)] <-
        oe[[pa("LOCF.count.step.", stepName)]][is.na(curColVec.old) & !is.na(curColVec.filled)] + 1 # increment LOCF count for all cells that changed from NA to not NA
    }
  }
  oe
}

# My reliability: SHR
ph.an.reliability.SHR <- function (vec.1, vec.2, ...) {
  vec.1 <- vec.1[!is.na(vec.1)]
  vec.2 <- vec.2[!is.na(vec.2)]
  sr.cor <- ph.an.cor(vec.1, vec.2, 
                      ...)
  ht <- an.sr.get.ht(sr.cor)
  stat.r.SBCorrected <- ph.an.formula.spearmanbrown.prophecy(ht$estimate, 2)
  sr.reliability <- anpriv.sr.create.reliability(
    dataNames=getSlot(sr.cor, "dataNames"),
    dataSource=getSlot(sr.cor, "dataSource"),
    desc=getSlot(sr.cor, "desc"),
    desc.set=getSlot(sr.cor, "desc.set"),
    ht=getSlot(sr.cor, "ht"),
    stat.r.SBCorrected=stat.r.SBCorrected,
    cronbach.alpha=ph.an.formula.cronbach.alpha(oe=as.data.frame(cbind(vec.1, vec.2)))
  )
  return (sr.reliability)
}

# Calculate split-half reliability, set
# Passed 2 oes with identical names
ph.an.reliability.SHR.set <- function (halfoe.1, halfoe.2, desc.set="reliability.SHR set", dataSource=as.character(match.call()[[2]]), ...) {
  srSet <- an.srSet.create(desc.set=desc.set)
  if (length(names(halfoe.1)) < 1)
    stop("ph.an.reliability.SHR(): too few columns in halfoe.1.")
  if (!identical(names(halfoe.1), names(halfoe.2)))
    stop("ph.an.reliability.SHR(): column names for halfoe.1 and halfoe.2 should be exactly the same. You could set names of one to equal the other, if they're just differing by")
  for (iColName in names(halfoe.1)[names(halfoe.1) != "sesType" & names(halfoe.1) != "p.id"]) {
    sr <- ph.an.reliability.SHR(halfoe.1[[iColName]],
                             halfoe.2[[iColName]],
                             dataNames=list(iColName, "half2"),
                             desc.set=desc.set,
                             dataSource=dataSource,
                             ...)
    srSet <- an.srSet.add(srSet, sr)
  }
  return (srSet)
}

# My t-test, independent samples
# Note - t.test defaults na.action to getOption('na.action'), which should be na.omit, leaving out any NAs from analyses.
ph.stat.t <- function (src.1, src.2, na.rm=T, ..., desc="", note="", srcNames=list(str.stripSpaces(deparse(match.call()[[2]])), str.stripSpaces(deparse(match.call()[[3]])))) {
  
  res <- t.test(src.1, src.2, ...)
  
  leveneRes <- ph.an.leveneTest.2vecs(src.1, src.2)
  p.levene <- leveneRes[["Pr(>F)"]][1]
  SDpool <- sqrt(((length(src.1)-1) * var(src.1, na.rm=na.rm) + (length(src.2)-1) * var(src.2, na.rm=na.rm)) /
                 (length(src.1) + length(src.2) - 2))

  rowoe <- data.frame(src.1=srcNames[[1]],
                      src.2=srcNames[[2]],
                      desc=desc,
                      sig=pa(ifelse(res$p.value < PH$STAT.P.TRIPLESTAR.CUTOFF, "***", 
                                    ifelse(res$p.value < PH$STAT.P.DOUBLESTAR.CUTOFF, "**",
                                           ifelse(res$p.value < PH$STAT.P.SINGLESTAR.CUTOFF, "*",
                                                  ifelse(res$p.value < PH$STAT.P.TREND.CUTOFF, "tnd", "")))),
                             ifelse(p.levene < PH$STAT.LEV.CUTOFF, "LEV", "")),
                      p=res$p.value[[1]],
                      M.diff=res$estimate[["mean of x"]] - res$estimate[["mean of y"]],
                      M.1=res$estimate[["mean of x"]],
                      M.2=res$estimate[["mean of y"]],
                      SDpool=SDpool,
                      d=(res$estimate[["mean of x"]] - res$estimate[["mean of y"]]) / SDpool,
                      CI.name.L=res$conf.int[1],
                      CI.name.U=res$conf.int[2],
                      levene=leveneRes[["Pr(>F)"]][1],
                      df=res$parameter[[1]],
                      stat.t=res$statistic[[1]],
                      note=note)
  CI.name.L.new <- pa("CI", attributes(res$conf.int)[["conf.level"]], "L")
  CI.name.L.new <- strs.replace(CI.name.L.new, "0.", "")
  names(rowoe)[names(rowoe)=="CI.name.L"] <- CI.name.L.new
  names(rowoe)[names(rowoe)=="CI.name.U"] <- strs.replace(CI.name.L.new, "L", "U")
  rowoe
}

# My one-sample t-test
ph.stat.t.1sample <- function (vec, mu=0, var.equal=F, ..., na.rm=T, desc="", desc.set="", dataSource="", dataNames=list(ph.str.stripSpaces(deparse(match.call()[[2]])))) {
  if (na.rm==T)
    vec <- vec[!is.na(vec)]
  res <- invisible(t.test(vec, mu=mu, var.equal=var.equal, ...))
  sr <- anpriv.sr.create.ttest.1sample(dataNames=dataNames,
                                       dataSource=dataSource,
                                       desc=desc,
                                       desc.set=desc.set,
                                       ht=list(res),
                                       mu=mu,
                                       sd=sd(vec))
  return (sr)
}

# Mean alternative: mean() with trimming
ph.an.mean.trimmed <- function (x, trim=.2, ...) {
  mean(x, tr=trim)
}

# Mean alternative: Just a wrapper for Wilcox's WRS winsorized mean
# Performance: Takes 10x as long as mean() without trimming
ph.an.mean.winsor <- function (x, trim=.2, ...) {
  WRS::win(x, tr=trim, ...)
}

# Wilcox-sourced robust effect size cid(), Cliff's method (Modern Statistics p.355)
# If you want a threat bias score, then x should be your neutral trials vector, y should be threat trials.
# (This matches the order of the subtraction for difference scores done as neutral - threat mean RTs)
ph.an.novelbias.ES.cliff <- function (x, y, ...) {
  cid(x=x, y=y, alpha=.05)$d # Note - cidv2 would be 15 times slower
}

# Cohen's d... this actually comes out quite different from bi.D. It's weird that the pooled SD thing makes such a difference
#   given that threat and neutral variances are overall the same. Maybe amongst people, they're meaningfully different, though I still don't see why it makes such a difference.
ph.an.novelbias.ES.cohen <- function (x, y, ...) {
  sd.pooled=sqrt(
    ( (length(x)-1) * var(x) + (length(y)-1) * var(y) ) /
      (length(x) + length(y) - 2))
  ( mean(x, ...) - mean(y, ...) ) /
    sd.pooled
}

# Wilcox-sourced robust effect size yuenv2() (Modern Statistics p.385)
# If you want a threat bias score, then x should be your neutral trials vector, y should be threat trials.
# This one is acting funny... it never goes up to 1.0 and actually gets lower if effect gets huge
ph.an.novelbias.ES.yuen <- function (x, y, trim=.2, nboot=100, ...) {
  yuenv2(x=x, y=y, tr=trim, alpha=.05, nboot=nboot, SEED=F)$Effect.Size
}

# Wilcox-sourced robust effect size akp.effect() (Modern Statistics p.385)
# If you want a threat bias score, then x should be your neutral trials vector, y should be threat trials.
# akp is just like cohen's d but with trimmed means and winsorized variance
ph.an.novelbias.ES.akp <- function (x, y, trim=.2, ...) {
  akp.effect(x=x, y=y, EQVAR=T, tr=trim)
}

# This is like akp, but it uses winsorizing instead of trimming on the mean
ph.an.novelbias.ES.akp.bothwin <- function (x, y, ...) {
  #[ base on source code of akp but use win mean instead of tmean ]
}

### You could do separate params for how much meantrimming and how much winsorizing is done on your variance... so you could try median but using winsorized var

# This code is mostly copied from Wilcox's akp, but it uses winsorizing instead of trimming on the mean
ph.an.novelbias.ES.bothwin <- function (x, y, trim=.2, EQVAR = TRUE, ...) {
  #[ base on source code of akp but use win mean instead of tmean ]
  tr <- trim
  x <- elimna(x)
  y <- elimna(y)
  n1 <- length(x)
  n2 <- length(y)
  s1sq = winvar(x, tr = tr)
  s2sq = winvar(y, tr = tr)
  winmean.x <- winmean(x, tr)
  winmean.y <- winmean(y, tr)
  
  spsq <- (n1 - 1) * s1sq + (n2 - 1) * s2sq
  sp <- sqrt(spsq/(n1 + n2 - 2))
  cterm = 1
  if (tr > 0) 
    cterm = area(dnormvar, qnorm(tr), qnorm(1 - tr)) + 2 * 
    (qnorm(tr)^2) * tr
  cterm = sqrt(cterm)
  if (EQVAR) 
    dval <- cterm * (winmean.x - winmean.y)/sp
  if (!EQVAR) {
    dval <- cterm * (winmean.x - winmean.y)/sqrt(s1sq)
    dval[2] = cterm * (winmean.x - winmean.y)/sqrt(s2sq)
  }
  dval
}

# Wilcox-sourced percentage bootstrap method, resorting to the p value as measure of ES.
# If you want a threat bias score, then x should be your neutral trials vector, y should be threat trials.
#   This is intended to be the order like traditionally "bias score = mean_neutral - mean_threat"
#   This ought to be revised... needs to be log or something to make p=.01 --> ES=.99 5 times as strong as .05
ph.an.novelbias.boot.pct.p <- function (x, y, trim=.2, nboot=599, ...) {
  res <- trimpb2(x=x, y=y, tr=trim, alpha=.05, nboot=nboot, SEED=F) # SEED=F is important!
  if (res$est.dif > 0)
    (1 - res$p.value)
  else
    -(1 - res$p.value)
}


ph.an.novelbias.bi.IAT <- function (x, y, ...) {
  ( mean(x, ...) - mean(y, ...) ) / sd(c(x, y))
}
ph.an.novelbias.bi.std <- function (x, y, ...) {
  ( mean(x, ...) - mean(y, ...) )
}

ph.an.novelbias.bi.std <- function (x, y, ...) {
  ( mean(x, ...) - mean(y, ...) )
}



# ------------------------------------ < Output and export functions > ------------------------------------

# Clears out the op/auto folder
ph.op.clrauto <- function () {
  unlink(pa('op/auto/', list.files('op/auto')))
}

# Outputs dataframe to csv for easy access
ph.op.oe <- function (oe, opTo="csv", fname="", run=T, runFile=run, cza=F, czsetNextRnd=T, na="") {
  if (opTo != "csv" & opTo != "tex")
    stop(pa("op(): Invalid opTo argument (\"", opTo, "\")"))
  wdSaved <- getwd()
  if (nchar(fname) == 0) {
    fname.core <- substitute(oe)
    fname.core <- paste(fname.core, collapse="\'")
    fname.parts <- c(path="op/auto",
                     core=ph.str.make.legal.filename(as.character(as.list(fname.core))),
                     ext=opTo)
    fname.parsed <- parse.fname.makeFromParts(fname.parts)
    if (!file.exists("op/auto")) {
      dir.create("op")
      dir.create("op/auto")
    }
    if (file.exists(fname.parsed[["whole"]])) {
      originalWhole <- fname.parsed[["whole"]]
      randomString <- ph.rndString(4)
      if (czsetNextRnd==T)
        assign("CZ.NEXT", randomString, envir=globalenv())
      fname.parsed[["core"]] <- pa(fname.parsed[["core"]], "-", randomString)
      fname.parsed <- parse.fname.makeFromParts(fname.parsed)
      ca("op(): Since \"", originalWhole, "\" exists already, creating new file with random suffix (", fname.parsed[["whole"]], ").")
    }
  } else {
    if (!ph.strs.detect(fname, opTo))
      fname <- pa(fname, ".", opTo)
    fname.parsed <- parse.fname(fname)
  }
  if (opTo=="csv") {
    outputFullPath <- fname.parsed[["whole"]]
    write.csv(oe, outputFullPath, row.names=F, na=na)
    PATH.OP.LAST <<- outputFullPath
    if (runFile==T)
      ph.run.file(outputFullPath)
  } else { # tex
    xt <- xtable(oe, type="latex")
    digits(xt) <- PH$FORMAT.DECIMALS
    ca(U.LATEX.BEGIN, print(xt), U.LATEX.END, file=fname.parsed[["whole"]])
    setwd(fname.parsed[["path"]])
    ph.run.shell(pa("\"", PATH.LATEX, "\\miktex\\bin\\pdflatex.exe\"", " ", fname.parsed[["proper"]]))
    outputFullPath <- str_replace(fname.parsed[["proper"]], "tex", "pdf")
    if (runFile==T)
      ph.run.file(outputFullPath)
    setwd(wdSaved)
  }
  if (cza==T)
    ph.cza(outputFullPath)
}
ph.op <- ph.op.oe


# Outputs table to csv with random filename in case op() can't make a valid filename out of it.
# Mainly for legacy.
# Pathprefix should specif a relative path plus the start of the filename to be used.
ph.op.oe.n <- function (oe, pathPrefix="opn/opn-", ...) {
  if (ph.isMode.skipOutput()) return (invisible(NULL))
  filename <- pa(pathPrefix, ph.rndString(4), ".csv")
  ca("opn(): Creating new file with random suffix (", filename, ").")
  write.csv(oe, filename, ...)
  ph.run.file(filename)
}
ph.opn <- ph.op.oe.n


# Outputs data frame to SPSS, creating a CSV and syntax based on your fname. This is awesome. 4/20/13
# You probably may want to run ph.oe.convert.numericToInteger() on your data frame first, to get rid of trailing .00 in SPSS.
# Behavior notes:
#  1. You must not use Excel to modify your CSV before opening via SPSS, because Excel messes it up by unquoting all text.
#  2. SPSS seems to have no way of understanding missing values in String variables, so any NAs from oe will become "" in SPSS
#     * That is, they will not properly become "MISSING" ... you should set the SPSS MISSING parameter in SPSS variable view to blank text ("" basically)
#  3. You should convert any date columns to proper R dates, with as.Date, before using this export
ph.op.oe.SPSS <- function (oe, fname="") {
  # Here's will need to turn dates into yyyy/mm/dd because the default output from R, of yyyy-mm-dd, SPSS won't take it
  # for classes, see list of basic types allowed in R data frames... also Date POSIXct...
  classes <- lapply(oe, class)
  classes <- lapply(classes, FUN=function (elem) { elem[[1]] }) # shorten to only first classes
  
  # Some classes translate into SPSS format with no changes to the data. These I call "easy" classes.
  ## Translate the easy class names into SPSS format names
  sTypes <- classes
  sTypes[classes == 'numeric']   <- 'F8.2'
  sTypes[classes == 'integer']   <- 'F8.0'
  sTypes[classes == 'logical']   <- 'F8.0'
  sTypes[classes == 'factor']    <- 'F8.0'
  
  toWrite.labels <- toWrite.nominals <- NULL
  # The "difficult" classes require changes to the actual data cols and specifying SPSS widths for Strings
  for (iCol in seq(ncol(oe))) {
    # Convert actual col data
    col.class <- classes[iCol]
    if (col.class == 'logical') {
      oe[[iCol]] <- as.integer(oe[[iCol]])
    } else if (col.class == 'POSIXct' | col.class == 'Date') {
      oe[[iCol]] <- format(oe[[iCol]], format="%Y/%m/%d")
    } else if (col.class == 'factor') {
      toWrite.labels <- c(toWrite.labels,
                          'VALUE LABELS '%_%names(oe)[iCol],
                          pa(1:length(levels(oe[[iCol]])), ' "', levels(oe[[iCol]]), '"'),
                          '.')
      toWrite.nominals <- c(toWrite.nominals, names(oe)[iCol])
      oe[[iCol]] <- as.integer(oe[[iCol]])
    }
    if (col.class == "character") {
      sTypes[iCol] <- 'A'%_%max(nchar(oe[[iCol]]))
    } else if (col.class == 'POSIXct' | col.class == 'Date') {
      sTypes[iCol] <- 'SDATE10'
    }
  }
  
  ph.op.oe(oe, fname=fname, run=F) # creates output file and sets PATH.OP.LAST to that path
  fconn <- file(PATH.OP.LAST%_%'.sps', 'w')
  PATH.OP.SPSS <- getwd()%_%'/'%_%PATH.OP.LAST
  PATH.OP.SPSS <- ph.strs.replace(PATH.OP.SPSS, '/', '\\')
  toWrite <- c(
    "*** This syntax file was generated by Phil's output to SPSS function (written 4/20/13) ***.",
    '******* Load from CSV file *******.',
    'GET DATA',
    '  /TYPE=TXT',
    '  /FILE="'%_%PATH.OP.SPSS%_%'"',
    '  /DELCASE=LINE',
    '  /DELIMITERS=","',
    "  /QUALIFIER='\"'",
    '  /ARRANGEMENT=DELIMITED',
    '  /FIRSTCASE=2',
    '  /IMPORTCASE=ALL',
    '  /VARIABLES=')
  
  
  toWrite <- c(toWrite,
               pa('    ', names(oe), ' ', sTypes),
               '.')
  
  toWrite <- c(toWrite,
               '******* Set value labels for factors *******.',
               toWrite.labels)
  
  toWrite <- c(toWrite,
               '******* Wrapping up *******.',
               'DATASET NAME DataSetFromR WINDOW=FRONT.',
               'VARIABLE WIDTH ALL (8).', 
               'VARIABLE LEVEL ALL (SCALE).',
               pa('VARIABLE LEVEL ', pa(collapse=' ', toWrite.nominals), ' (NOMINAL).'))
  
  writeLines(toWrite, con=fconn)
  close(fconn)
}

# opm... opm is like output to my stream, which will save an oe of a message, which will often be just 1 string (pa'd by caller)
#  it will be able to output to the console (mmoc) or to an opn excel (mmoe)
#  param t means show timestamp, param s means silent
# The ... is for extra print args like "row.names=F, digits=4"... they won't have any effect on the opm that is stored; they'll only affect the opm that's shown to user when s=F

opm <- function (item, t=F, s=F, timer=F, row.names=F, tight=F, margin.top.oe=(tight+1)%%2, margin.bot.oe=(tight+1)%%2, ...) { 
  if (ph.isMode.skipOutput()) return (invisible(NULL))
  if (missing(item)) {
    ca(ph.timeStamp(time=T, timer=timer))
    return (invisible())
  }
  if (!exists("opm.msgList"))
    opmc()
  msgList <- get("opm.msgList", envir=baseenv())
  assign("opm.msgList", c(msgList, list(list(timeStamp=ph.timeStamp(time=T, timer), item=item, tight=tight))), envir=baseenv())
  if (s == FALSE) {
    if (t==TRUE | timer==T)
      ca(ph.timeStamp(time=t, timer=timer), " ")
    if (length(item) > 1 | is.list(item)) {
      #       ca(":\n\n")
      if (is.data.frame(item)) {
        ca(rep('\n', margin.top.oe))
        print(item, row.names=row.names, ...)
        ca(rep('\n', margin.bot.oe))
      } else {
        ca(rep('\n', margin.top.oe))
        print(item, ...)
        ca(rep('\n', margin.bot.oe))
      }
    } else {
      #       ca(rep('\n', 1))
      cal(item)
      #       ca(rep('\n', 1))
    }
  }
}
ph.opm <- opm

# Clears opm / creates a new empty stream
opmc <- function () {
  assign("opm.msgList", list(), envir=baseenv())
  ca("opmc(): Created new opm stream.\n")
  assign("opm.lastopmsTime", Sys.time() - 1, envir=baseenv())
}
ph.opmc <- opmc


# Outputs last n messages
opms <- function (n=-1, to="console", open=T) {
  if (ph.isMode.skipOutput()) return (invisible(NULL))
  recentOnly <- ifelse(n == -1, TRUE, FALSE)
  nOut <- 0
  if (!exists("opm.msgList") | !exists("opm.lastopmsTime"))
    return ("opms(): No messages to output.")
  switch (to,
          "console" = filename <- "console",
          "txt" = filename <- "opms.txt",
          "csv" = filename <- "opms.csv",
          filename <- to
  )
  if (substr(filename, str_length(filename)-2, str_length(filename)) == "csv") {
    csvMode <- TRUE
  } else {
    csvMode <- FALSE
  }
  if (filename != "console") {
    unlink(filename)
    if (!file.exists(file.path("op"))) {
      dir.create("op")
    }
    sink(file.path("op", filename))
  }
  msgList <- get("opm.msgList", envir=baseenv())
  if (length(msgList) == 0)
    return ("opms(): No messages to output.")
  if (n > length(msgList) | recentOnly)
    n <- length(msgList)
  #   if (substr(filename, str_length(filename)-2, str_length(filename)) == "csv")
  #     ca("Time, Elapsed, Unit, Item\n")
  for (i in 1:n) {
    msg <- msgList[[length(msgList) - n + i]]
    if (recentOnly & strptime(ph.strs.extractBetween(msg$timeStamp, "(", "|"), format="%H:%M:%S") < get("opm.lastopmsTime", envir=baseenv()))
      next
    nOut <- nOut + 1
    item <- msg$item
    if (csvMode==T) {
      cat.for.csv(item=item, tight=msg$tight)
    } else {
      if (is.list(item) | length(item) > 1) {
        if (msg$tight==F) cal()
        print(item)
        if (msg$tight==F) cal()
      } else {
        cal(item)
      }
    }
  }
  toCat <- pa("opms(): Outputted ", nOut, " items (", n, " available)")
  if (csvMode==T)
    toCat <- ph.up.cat.str.csv.prep(toCat)
  ca(toCat)
  
  if (filename != "console") {
    sink()
    if (open == T)
      ph.run.file(file.path("op", filename))
  }
  assign("opm.lastopmsTime", Sys.time(), envir=baseenv())
}
ph.opms <- opms

# opm convenience aliases
opmsa <- function (...) opms(n=999, ...)
ph.opmsa <- opmsa


# ------------------------------------ < All other functions > ------------------------------------

ph.isMode.skipOutput <- function() return (PH$DISABLE.SOME.OUTPUT)

# From function curfnfinder() in post http://stackoverflow.com/questions/7307987/logging-current-function-name
ph..curFunctionName <- function (skipframes=0, skipnames="(FUN)|(.+apply)|(replicate)|.fun",
                              retIfNone="[Not currently in a function]", retStack=FALSE, extraPrefPerLevel="\t")
{
  prefix<-sapply(3 + skipframes+1:sys.nframe(), function(i){
    currv<-sys.call(sys.parent(n=i))[[1]]
    return(currv)
  })
  prefix[grep(skipnames, prefix)] <- NULL
  prefix<-gsub("function \\(.*", "do.call", prefix)
  if(length(prefix)==0)
  {
    return(retIfNone)
  }
  else if(retStack)
  {
    return(paste(rev(prefix), collapse = "|"))
  }
  else
  {
    retval<-as.character(unlist(prefix[1]))
    if(length(prefix) > 1)
    {
      retval<-paste(paste(rep(extraPrefPerLevel, length(prefix) - 1), collapse=""), retval, sep="")
    }
    return(retval)
  }
}

# From http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
ph..ls.objects <- function (pos = 1, pattern, order.by, decreasing=FALSE, head=FALSE, n=5) {
  if (ph.isMode.skipOutput()) return (invisible(NULL))
  napply<- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(print(object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# ph..ls.objects() shortcut
ph.lsos <- function(..., n=10) ph..ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)

# Opens in Windows Explorer all folders that you have a "PATH." variable for, such as PATH.HELPER and PATH.USER.DATA
ph..paths <- function () {
  if (ph.isMode.skipOutput()) return (invisible(NULL))
  pathVarVec <- ls(envir=globalenv())[which(ph.strs.detect.vec(substr(ls(envir=globalenv()), 1, 5), "PATH."))]
  pathVarVec <- pathVarVec[pathVarVec != "PATH.LATEX" & pathVarVec != "PATH.DBOX"]
  pathVec <- vector()
  for (iPathVar in pathVarVec)
    pathVec <- c(pathVec, get(iPathVar))
  pathVec <- union(pathVec, getwd())
  for (iPath in pathVec) {
    ph.run.file(iPath)
  }
}

# Shortcut for as.factor(as.character(vec))
ph.as.fac <- function (vec) as.factor(as.character(vec))

# Asserts that there are no NA's in passed vector or data frame... errors if there are
ph.assertNoNAs <- function (data) {
  if (is.vector(data) | is.data.frame(data)) {
    if (is.data.frame(data) & nrow(data) == 0)
      return (invisible(NULL))
    if (ph.freqTRUE(is.na(data)) != 0)
      stop("ph.assertNoNAs(): You do have NAs in: \"", substitute(data), "\"! Halting right now. -Phil")
  } else {
    stop("assertNoNAs: I'm not sure how to handle what you passed me--it's not a vector or data frame. -Phil")
  }
}

# My benchmarking function, makes random numbers and rbinds in a big dataframe, base R only
# Results are 17.8s on my Home i7, 15.2 on my new laptop i7 plugged in, and 34.5s when it's unplugged
ph.benchmark <- function (nReps=20, nRows=50, nCols=500) {
  cat("'user' shows time elapsed (seconds) attributable to executing our benchmarking code\n")
  print(system.time(
{
  for(iRep in seq(nReps)) {
    x <- rnorm(nCols)
    oe <- as.data.frame(t(x))
    for (i in seq(nRows)) {
      x <- rnorm(nCols)
      oe <- rbind(oe, t(x))
    }}}))
}

# ca() shortcuts
ca <- function(..., sep="") {
  if (ph.isMode.skipOutput()) return (invisible(NULL))
  cat(..., sep=sep)
}
ph.ca <- ca
cal <- function(...) ca(..., "\n")
ph.cal <- cal

# cat a data frame as comma separated
ph.up.cat.oe.csv <- function(oe, digits=15, scipen=4) {
  if (nrow(oe) == 0 | ncol(oe) == 0) {
    return (print(oe))
  }
  cal('"', pa(collapse='","', names(oe)), '"')
  for (iRow in 1:nrow(oe)) {
    for (iCol in 1:ncol(oe)) {
      toOut <- oe[[iRow, iCol]]
      switch(class(toOut),
             "factor" = {
               toOut <- as.character(toOut)
             },
             "numeric" = {
               toOut <- format(toOut, scientific=scipen, digits=digits)
             },
             "character" = {
               toOut <- '"'%_%toOut%_%'"'
             })
      ca(toOut, ",")
    }
    cal()
  }
}

# Prepares string for csv'ing by putting quotes around it and any line breaks within it
ph.up.cat.str.csv.prep <- function (text) {
  text <- ph.strs.replace(text, '\n', '"\n"')
  text <- '"'%_%text%_%'"'
  text
}

# Cats for a csv... it's assuming you're sinking into a csv.
cat.for.csv <- function (item, tight=F) {
  if (is.data.frame(item) == TRUE) {
    if (tight==F) cal()
    ph.up.cat.oe.csv(item)
    if (tight==F) cal()
  } else {
    if (is.list(item) | length(item) > 1) {
      if (tight==F) cal()
      print(item)
      if (tight==F) cal()
    } else {
      if (is.character(item))
        item <- ph.up.cat.str.csv.prep(item)
      cal(item)
    }
  }
}

# Assigns to cz, clears whole thing and assigns
ph.czass <- function (newcz) {
  assign("CZ", newcz, envir=globalenv())
}

# Gets CZ from global env
ph.czget <- function () {
  get("CZ", envir=globalenv())
}

# Creates a "codezip" (my idea) as an archive of all code that I need for a specific analyses, even a tiny variant of analysis... plus sunk output files too hopefully
ph.cz <- function (files=ph.czget(), zipName=NULL, folder=pa(PATH.DBOX, "/R/codezips")) {
  #   files <- c(".Rhistory", files)
  
  files <- sapply(files, function (fname) { return (pa("\"", fname, "\"")) })  
  if (missing(zipName)) {
    go <- T
    while (go==T) {
      if (exists("CZ.NEXT", envir=globalenv())) {
        zipName <- get("CZ.NEXT", envir=globalenv())
        rm("CZ.NEXT", envir=globalenv())
      } else {
        zipName <- ph.rndString(4) # 1.7 million possibilities
      }
      if (!file.exists(file.path(folder, zipName)))
        go <- F
    }
  } else {
    if (file.exists(file.path(folder, zipName)))
      file.remove(file.path(folder, zipName))
  }
  shell(pa(PATH.DBOX,
           "/Utils/7Zip/App/7-Zip/7z.exe a -tzip -mx9 ",
           file.path(folder, zipName),
           " ",
           paste(files, collapse=" ")))
  ca("czipped into: '", zipName, "'")
  assign("CZ.LAST", file.path(folder, zipName), envir=globalenv())
}

# Adds new file path to CZ *if* it's not already in there
ph.cza <- function (newFile) {
  files <- ph.czget()
  if (!any(files == newFile)) {
    files <- c(files, newFile)
    assign("CZ", files, envir=globalenv())
  }
}

# Replaces most recent addition to czip
ph.czrepl <- function(newFile) {
  files <- ph.czget()
  files[length(files)] <- newFile
  assign("CZ", files, envir=globalenv())
}

# Runs last filepath in CZ
ph.czrunl <- function () {
  ph.run.file(ph.czget()[[length(ph.czget())]])
}

# Open czip's folder (path)
ph.czp <- function() {
  ph.run.file(pa(PATH.DBOX, "/R/codezips"))
}

# Removes most recent czip, or the specified one
ph.czrm <- function (toRemove=NULL) {
  if (is.null(toRemove)) {
    toRemove <- CZ.LAST
  } else {
    toRemove <- pa(PATH.DBOX, "/R/codezips/", toRemove)
  }
  ret <- file.remove(pa(toRemove, ".zip"))
  if (ret == T)
    ca("Successfully removed: ", pa(toRemove, ".zip"))
}

# Unzips to a folder of that name, still in the codezips folder... you can just manually delete all folders once in a while from there
ph.czu <- function(zipName, folder=pa(PATH.DBOX, "/R/codezips"), destPath=pa(folder, "\\", zipName)) {
  zPath <- pa(PATH.DBOX, "/Utils/7Zip/App/7-Zip/7z.exe")
  zPath <- str_replace_all(zPath, "/", "\\\\")
  folder <- str_replace_all(folder, "/", "\\\\")
  destPath <- str_replace_all(destPath, "/", "\\\\")
  zipName <- str_replace_all(zipName, "/", "\\\\")
  if (!file.exists(file.path(folder, pa(zipName, ".zip")))) {
    ca("unczip(): Couldn't find your '", zipName, ".zip'")
  } else {
    shell(pa(zPath, " x -y -o",
             destPath,
             " ",
             folder, "\\", zipName, ".zip"))
  }
  ph.run.file(destPath)
}

# Creates blank data frame of however many columns
ph.data.frame.empty <- function (colNames, nrow=0) {
  oe <- data.frame(matrix(nrow=nrow, ncol=length(colNames)))
  names(oe) <- colNames
  oe
}

# My descriptives
ph.dsc <- function (x, hist=F, show.NA, ...) {
  #   ca("mean ", mean(x, ...), "; median ", median(x, ...), "; sd ", sd(x, ...))
  ph.outoe.make(c("n.total", "n.NA", "mean", "median", "sd", "min", "max" ))
  n.total <- length(x)
  n.NA <- sum(is.na(x))
  x <- x[!is.na(x)]
  ph.outoe.push(data.frame(n=n.total,
                        n.NA=n.NA,
                        mean=mean(x, ...),
                        median=median(x, ...),
                        sd=sd(x, ...),
                        min=min(x, ...),
                        max=max(x, ...)))
  if (hist==T)
    hist(x, breaks=50)
  ph.outoe.get()
}

# Reads all CSVs in a folder into a data frame
ph.file.stackCSVs <- function (src.path,
                               src.filenames.match='.csv', # can be any string, e.g., '_subject-', and we'll only grab filenames containing this
                               print.filenames=T,
                               addCol.src.folder=F,   # creates col in the oe containing the folder path of the file that the row came from
                               addCol.src.filename=F, # ... same for file name
                               addCol.rowNum=F,
                               header=T,
                               stringsAsFactors=T,
                               ... # passed to read.csv (which calls read.table with ... too)
                               ) {
  oeList <- NULL
  if(print.filenames)
    opm('ph.file.stackCSVs(): Reading CSVs from "' %_% src.path %_% '":')
  for(nextFile in list.files(full.names=T, src.path)) {
    filename.proper <- strs.extractAfter(nextFile, '/', fromLast=T)
    if (strs.detect(filename.proper, src.filenames.match)) {
      if(print.filenames)
        opm(filename.proper)
      nextFileoe <- read.csv(nextFile, header=header, stringsAsFactors=stringsAsFactors)
      if (addCol.src.folder)
        nextFileoe$src.folder <- strs.replace(nextFile, filename.proper, '')
      if (addCol.src.filename)
        nextFileoe$src.filename <- filename.proper
      if (addCol.rowNum)
        nextFileoe$rowNum <- seq(1, nrow(nextFileoe))
      oeList <- c(oeList, list(nextFileoe))
    }
  } 
  oe <- rbind.fill(oeList)
  oe
}

# Shortcuts for formatting numerics.
ph.fmat.sci <- function (x) {
  format(x, scientific=T)
}
ph.fmat.dig <- function (x) {
  format(x, digits=22) # 22 is maximum
}
ph.fmat.scidig <- function (x) {
  format(x, digits=22, scientific=T)
}

# Wrapper for prettyR's freq function, to enable the capability of global output disabling
ph.freq <- function (...) {
  if (ph.isMode.skipOutput())
    return (invisible(NULL))
  else
    get("freq", envir=getNamespace("prettyR"), inherits=F)(...)
}

# Display frequencies, setting a max length on levels, in case the names are too long
ph.freq.short <- function(var, len=11) {
  ph.freq(substr(var, 1, len))
}

# Return frequency of needle in vector haystack. Warn if any NAs exist.
# (Btw, you used as.character(needle) to make sure)
ph.freqOf <- function (haystack, needle) {
  if (sum(is.na(haystack)) > 0)
    ca("ph.freqOf(): Note - ", as.character(substitute(haystack)), " has NAs: ", sum(is.na(haystack)), "\n")
  length(haystack[haystack==needle])
}
# Convenience methods for freqOf
ph.freqTRUE <- function (haystack) {
  ph.freqOf(haystack, TRUE)
}
ph.freqFALSE <- function (haystack) {
  ph.freqOf(haystack, FALSE)
}

# Return a logical vector with all duplicate entries INCLUDING all originals
ph.getAllDupes <- function (vec) {
  duplicated(vec) | duplicated(vec, fromLast=T)
}



# Takes a list with 2 (or more) levels
# Say you have a list of 5 lists, each has element A and B
# Use this function to create a list of 5 As (or 5 Bs, just specify which elementName you want, like "A" or "B")
ph.list.list.makeListFromLevel2Element <- function (list.list, elementIndex=NA, elementName=NA) {
  if ( (!is.na(elementIndex) & !is.na(elementName)) |
         (is.na(elementIndex) & is.na(elementName)) )
    stop("list.list.makeListFromLevel(): You must set either elementIndex OR elementName (not both, not neither).")
  retList <- vector("list", length(list.list))
  for (i in 1:length(list.list)) {
    listItem.level1 <- list.list[[i]]
    subscript <- ifelse(is.na(elementIndex), elementName, elementIndex)
    retList[[i]] <- listItem.level1[[subscript]]
  }
  retList
}


# Output all the data into csvs, so you can look at the tree in windows explorer
# Can be applied to 1-level lists as well as multi-level
ph.list.list.op.tree <- function (list.list, opTo="csv", folder=".") {
  if (ph.strs.detect.suffix(folder, "\\")) # Remove ending slash because file.exists needs that removed
    folder <- substr(folder, 1, nchar(folder)-1)
  if (!file.exists(folder))
    dir.create(folder)
  for (iElem in 1:length(list.list)) {
    elem <- list.list[[iElem]]
    if (class(elem)=="list") {
      newFolderName <- names(list.list)[iElem]
      if (newFolderName == "")
        newFolderName <- pa("elem.", iElem)
      ph.list.list.op.tree(list.list=elem, opTo=opTo, folder=pa(folder, "\\", newFolderName))
    }
    else {
      elemName <- names(list.list)[iElem]
      if (elemName == "")
        elemName <- pa("elem.", iElem)
      lastBreadcrumb <- substr(folder, ph.strs.pos(folder, "\\")+1, nchar(folder))
      cal("Creating file: ", pa(folder, "\\", lastBreadcrumb, "-", elemName, ".", opTo))
      op(getOutput(elem),
         opTo=opTo,
         fname=pa(folder,
                  "\\",
                  ph.str.make.legal.filename(lastBreadcrumb),
                  "-",
                  elemName),
         runFile=F)
    }
  }
}

# Throughout a whole list tree, eliminates any lists of length 1 by promoting the 1 element to replace its parent list
ph.list.list.promoteSingleElementLists <- function (lst) {
  if (!is.list(lst) | is.data.frame(lst))
    return (lst)
  for (iElem in 1:length(lst)) {
    if (is.list(lst[[iElem]])) {
      if (!is.data.frame(iElem)) {
        if (length(lst[[iElem]]) == 1)
          lst[[iElem]] <- lst[[iElem]][[1]]
        else
          lst[[iElem]] <- ph.list.list.promoteSingleElementLists(lst[[iElem]])
      }
    }
  }
  lst
}

# Reorients a list of (lists with indexable elements)
# For example, if you start with 10 lists each with 3-length lists, this will return 3 lists each with 10-length lists
ph.list.list.reorient <- function (list.list) {
  jLength <- length(list.list[[1]])
  retList.list <- vector("list", jLength)
  for (j in 1:jLength) {
    retList.list[[j]] <- list()
  }
  for (iRow in 1:length(list.list)) {
    for (jRow in 1:jLength) {
      retList.list[[jRow]] <- c(retList.list[[jRow]], list.list[[iRow]][[jRow]])
    }
  }
  retList.list
  # test code: ph.list.list.reorient(list(list(list(2), list(20)), list(list(3), list(30))))
}



# Takes a list of oes--or a list of anything whose getOutput() returns an oe--and returns one bound-together oe.
# Right after doing getOutput on a given element to get an oe, the subset expression is performed.
# So, in your subsetxp, use the context (column/variable names) of the source oe.
# Put your subsetxp in quotes, and escape your quotes when needed.
ph.list.oe.makeoeFromSubset <- function(list.oe, subsetxp) {
  return (do.call(rbind, lapply(list.oe, FUN=function (item) {
    output <- getOutput(item)
    if (is.data.frame(output))
      return (subset(output, eval(parse(text=subsetxp))))
  }
  )))
}

# rbinds together a list of oes
ph.list.oe.makeoe <- function(list.oe) {
  ph.list.oe.makeoeFromSubset(list.oe, TRUE)
}

ph.list.oe.matchColPrintWidths <- function(list.oe, digits=getOption("digits"), nsmall=0) {
  # Private function: Returns nchar of a format 
  nchar.format<- function (x, digits=getOption("digits"), nsmall=0) {
    return (nchar(format(x, digits=digits, nsmall=nsmall)))
  }
  
  l <- list()
  l[[1]] <- soe
  l[[2]] <- soe
  l
  
  nColsVec <- sapply(l, ncol)
  if (!all(nColsVec==nColsVec[1]))
    stop("ph.list.oe.matchColPrintWidths(): Your oes have differing nCols.")
  
  for (i in seq_along(list.oe)) {
    oe <- list.oe[[i]]
    for (j in seq_along(oe)) {
      targetLength <- max(nchar(
        as.character(
          ifelse(is.numeric(oe[[j]]), format(oe[[j]], digits=digits, nsmall=nsmall), oe[[j]])
        )))
      names(oe)[j] <- ph.str.padToLength(names(oe)[j], targetLength=targetLength)
    }
    for (j in seq_along(oe)) {
      #       # NEXT ACTION: Here, compare to lengths of all other oes' columns, by number. Add initial check of nCols, that should match for all oes.
      #       max(lapply(list.oe, FUN=nchar.format,
      #                           digits=digits,
      #                           nsmall=nsmall
      #                  ))
    }
    list.oe[[i]] <- oe
  }
  list.oe
}


# From a list of protos, returns a vector of returns from whatever method
ph.list.proto.methodCall.vecReturn <- function (protoList, methodName) {
  retVec <- vector()
  for (iProto in protoList) {
    iProto[["methodName"]] <- methodName
    retVec <- c(retVec, with(iProto, get(methodName))(iProto))
    iProto[["methodName"]] <- NULL
  }
  retVec
}


# Removes elements from a list...
# This will not remove any unkeyed elements (or elements with key ""); it only affects named elements
ph.list.remove.elements <- function (lst, ..., silent=T) {
  keysToKeep <- c("", ph.vec.remove.elements(names(lst), ..., silent=T))
  keysToRemove <- ph.listAbsences(keysToKeep, names(lst))
  for (iKey in keysToRemove)
    lst[iKey] <- NULL
  if (silent==F)
    ca("ph.list.remove.elements(): Removed ", nColsBefore - length(keysToKeep), " elements, ", length(unique(colNumsToKeep)), " remaining\n")
  lst
}

# List needles that aren't in haystack
ph.listAbsences <- function (haystack, needles) {
  absenceVec <- NULL
  for (needle in needles) {
    if (!(needle %in% haystack))
      absenceVec <- c(absenceVec, needle)
  }
  absenceVec
}

# List elements of first vector that ARE ALSO in second vector
ph.listPresences <- function (haystack, needles) {
  presenceVec <- NULL
  for (needle in needles) {
    if (needle %in% haystack)
      presenceVec <- c(presenceVec, needle)
  }
  presenceVec
}

# Does a load, into globalenv(), just with more info
ph.load.withmsg <- function(file, envir=globalenv()) {
  loadedNames <- load(file, envir=envir)
  opm("Loaded '"%_%file%_%"' containing '"%_%pa(loadedNames, collapse=',')%_%"' (wd: '"%_%getwd()%_%"')")
}

# mfm = myFormat... surround a string with suitable header marks
ph.mfm <- function (type, ..., margin.top=0, margin.bottom=0, flankString="|") {
  switch (type,
          "h1" = {
            nTailStars <- 120 - length("---------------- ") - nchar(pa(...)) - 1
            toret <- pa("---------------- ", ..., " ", paste(rep("-", nTailStars), collapse=""))
          },
          "h1.5" = {
            toret <- pa("_-_-_-_-_-_-_-_- ", ..., " -_-_-_-_-_-_-_-_")
          },
          "h2" = {
            toret <- pa("___________ ", ..., " ___________")
          },
          "h2.5" = {
            toret <- pa("*_*_*_*_* ", ..., " *_*_*_*_*")
          },
          "h3" = {
            toret <- pa("****** ", ..., " ******")
          },
          "h4" = {
            toret <- pa("___ ", ..., " ___")
          },
          "h1note" = {
            toret <- pa("********** ", ...)
          },
          "h2note" = {
            toret <- pa("******** ", ...)
          },
          "h3note" = {
            toret <- pa("****** ", ...)
          },
          "h4note" = {
            toret <- pa("*** ", ...)
          },
          "flank" = {
            toret <- pa(flankString, ..., flankString)
          },
          "flank right" = {
            toret <- pa(..., flankString)
          },
          "flank left" = {
            toret <- pa(flankString, )
          },
{
  stop("ph.mfm(): You gave me an invalid formatting type: \"", type, "\"")
}
  )
  
  if (margin.top > 0)
    toret <- pa(pa(rep("\n", margin.top), collapse=""), toret)
  if (margin.bottom > 0)
    toret <- pa(toret, pa(rep("\n", margin.bottom), collapse=""))
  toret
}
ph.mfm.vcz <- Vectorize(ph.mfm)


# Display names(oe) but quoted and separated by commas, useful sometimes for copying and pasting
ph.names.cat <- function (oe) {
  theList <- ''
  for (iVar in names(oe)) {
    theList <- paste(theList, '"', iVar, '", ', sep='')
  }
  ca(theList)
}



# Renames oe's columns, using the syntax of a recode from the 'car' package
# @returns list with $oe, $records, and $colLabels
ph.oe.cols.rename <- function (oe, recodes, colLabels=NULL, silent=F) {
  recodeRet <- ph.vec.recode(vec=names(oe), recodes=recodes, valueLabels=colLabels, silent=silent)
  names(oe) <- recodeRet$vec
  return (list(oe = oe, records = recodeRet$records, colLabels = recodeRet$valueLabels))
}

# TO DO: might want to rewrite this with reshape2::melt() and cast(), http://goo.gl/EqLgv ... I dunno, maybe this one's good already
## Converts from long (like sesoe) to wide (like poe) format
# Should work for any situation of multiple observations per subject, not just sessTypes / timepoints. Made 4/18/13
# Param: cn.sessType, column name for a col that contains timepoints or sesstypes. Will not exist in output oe.
# Param: cns.P specifying col names that are persistent across all sessTypes (like P signup data). These must be identical in all sessType entries
#   cn.p.id is automatically included. It should error if you forget to include cns.P cols that aren't persistent
ph.oe.convert.longToWide <- function (oe, cn.p.id, cns.P=NULL, cn.sessType, sessTypes='[ALL]', separator='.') {
  if (!all(cns.P %in% names(oe)))
    stop("ph.oe.convert.longToWide(): All cns.P should be column names in your oe... I didn't find all of them in there.")
  stopifnot(cn.p.id %in% names(oe) & cn.sessType %in% names(oe))
  cns.P <- union(cns.P, cn.p.id)
  cns.measures <- setdiff(names(oe), cns.P) # "measures" should be any cols not persistent throughout sessTypes
  cns.measures <- cns.measures[cns.measures != cn.sessType]
  
  oe[[cn.sessType]] <- as.character(oe[[cn.sessType]])
  if (sessTypes[1] == '[ALL]') {
    sessTypes <- unique(sort(oe[[cn.sessType]]))
  } else {
    sessTypes <- as.character(sessTypes)
    if (!all(sessTypes %in% unique(oe[[cn.sessType]])))
      stop("ph.oe.convert.longToWide(): Your sessTypes should all appear in your sessType column '"%_%cn.sessType%_%"' but they don't")
    stopifnot(length(sessTypes) == length(unique(oe[[cn.sessType]])))
  }
  
  poe <- oe[cns.P]
  poe <- unique(poe)
  if (length(unique(poe[[cn.p.id]])) < nrow(poe))
    stop("ph.oe.convert.longToWide(): Looks like you have some data in the cns.P cols that aren't identical across all those Ps' rows")
  
  for (i in seq_along(sessTypes)) {
    # Make the toe with just p.id and the measures, appending sessType name to colnames
    toe <- oe[oe[[cn.sessType]] == sessTypes[i],
              names(oe) == cn.p.id | names(oe) %in% cns.measures]
    if (any(duplicated(toe[[cn.p.id]]) == TRUE))
      stop("ph.oe.convert.longToWide(): There's more than one row for sessType "%_%sessTypes[i]%_%" for some P")
    nms <- names(toe)
    nms[nms %in% cns.measures] <- paste(sep="", nms[nms %in% cns.measures], separator, sessTypes[i])
    names(toe) <- nms
    poe <- merge(poe, toe, by=cn.p.id, all.x=T, all.y=F)
  }
  as.data.frame(poe)
}

# Converts wide to long format, intended for time and repeated measures
# It will convert to long format any cols with any ext as a suffix
#   Any cols that don't have an ext as a suffix are treated same as 'p.' cols: They'll be same value on every row of a given P
# Note - it's ok if not all exts you pass are actually present in colnames
# melt() is from reshape package... 
# Cool reshape guide (see the diagram a ways down especially): http://www.r-statistics.com/2012/01/aggregation-and-restructuring-data-from-r-in-action/
ph.oe.convert.wideToLong <- function (oe,
                                   exts,
                                   timeColName = 'time',
                                   cns.keepFirst = c('p.id', 'p.condition', timeColName), # keeps first IF PRESENT, so you can put more in arg
                                   convertNumericOnly = T) {
  cns.exted <- ph.oe.ext.getExtedColnames(oe, exts)
  cns.persistent <- names(oe)[(!(names(oe) %in% cns.exted))]
  if (convertNumericOnly==T)
    cns.exted <- cns.exted[lapply(oe[cns.exted], is.numeric) == TRUE]
  oe <- oe[c(cns.persistent, cns.exted)]
  library(reshape)
  moe <- melt(oe, id.vars=cns.persistent)
  moe$variable <- as.character(moe$variable)
  moe[[timeColName]] <- NA
  for (ext in exts) {
    extRowLogi <- ph.strs.detect.suffix.vec(moe$variable, ext)==T
    if (sum(extRowLogi) == 0)
      next
    moe[extRowLogi, timeColName] <- ext
    moe[extRowLogi, 'variable'] <- strs.replace.suffix.vec(moe[extRowLogi, 'variable'], ext, '')
  }
  coe <- cast(moe, ...~variable)
  nms <- names(coe)
  cns.keepFirst <- intersect(cns.keepFirst, nms)
  cns.wide <- nms[!(nms %in% cns.keepFirst) & !(nms %in% cns.persistent)]
  cns.persistent.notFirst <- cns.persistent[!(nms %in% cns.keepFirst)]
  coe <- coe[c(cns.keepFirst, sort(cns.wide), sort(cns.persistent.notFirst))]
  coe[[timeColName]] <- as.factor(coe[[timeColName]])
  as.data.frame(coe)
}


# Make (or update) difference scores
# Or give exts and stems, like exts=c("w0", "w4"), stems=c("LSAS", "PSWQ")
# For suffix, you may want like ".w0_w4.D"
# The D score is calculated as post minus pre.
ph.oe.Ds.make <- function (oe, stems='[ALL]', ext.pre, ext.post, suffix='.'%_%ext.pre%_%'_'%_%ext.post%_%'.D', silent=F) {
  if (stems[1] == '[ALL]') {
    validStems <- character()
    for (nm in names(oe)) {
      if (ph.strs.detect.suffix(nm, ext.pre)==T) {
        stem.found <- ph.strs.replace(nm, ext.pre, "")
        if (stem.found%_%ext.post %in% names(oe)) {
          if (is.numeric(oe[[stem.found%_%ext.pre]]) & 
                is.numeric(oe[[stem.found%_%ext.post]])) {
            validStems <- c(validStems, stem.found)
          }
        }
      }
    }
    stems <- validStems
    opm("Found stems: "%_%pa(stems, collapse=', '))
  }
  for (stem in stems) {
    if (!(stem%_%ext.pre %in% names(oe) & stem%_%ext.post %in% names(oe))) {
      stop("ph.oe.Ds.make(): stem '"%_%stem%_%"' doesn't exist in oe with both exts.")
    }
    if (!is.numeric(oe[[stem%_%ext.pre]]) | 
          !is.numeric(oe[[stem%_%ext.post]])) {
      stop("ph.oe.Ds.make(): stem '"%_%stem%_%"' is invalid, as its cols aren't numeric")
    }
    oe[[stem%_%suffix]] <- oe[[stem%_%ext.post]] - oe[[stem%_%ext.pre]]
  }
  oe
}


# Extends a stem variable name, returning any valid column names
ph.oe.ext <- function (oe, stems, exts, onlyReturnExtsThatAllStemsHave=F, returnOnlyValidColNames=T) {
  if (onlyReturnExtsThatAllStemsHave==T)
    validExts <- ph.oe.ext.validExts(oe=oe, stems=stems, exts=exts)
  validColNames <- character()
  for (stem in stems) {
    if (onlyReturnExtsThatAllStemsHave==F)
      validExts <- ph.oe.ext.validExts(oe=oe, stems=stem, exts=exts)
    if (!is.null(validExts))
      validColNames <- c(validColNames, stem%_%validExts)
  }
  validColNames
}

# Checks which of passed exts exist for the stems, returning any exts that exist for all of them
ph.oe.ext.validExts <- function (oe, stems, exts) {
  validExts <- character()
  for (iExt in exts) {
    if (stems[1]%_%iExt %in% names(oe))
      validExts <- c(validExts, iExt)
  }
  invalidatedExts <- character()
  if (length(stems) > 1) {
    for (stem in stems[2:length(stems)]) {
      for (ext in validExts) {
        if (!(stem%_%ext %in% names(oe)))
          invalidatedExts <- union(invalidatedExts, ext)
      }
    }
  }
  validExts <- validExts[!(validExts %in% invalidatedExts)]
  validExts
}

# Does sapply on columns of an oe made from ext'ing the stem
# Tries to return an oe in long format, with the exts as first column
# With na.listwise on, it's like sapply(na.omit(oe[ph.oe.ext(oe)]), fun)... might be better to just use that na.omit manually
ph.oe.ext.apply <- function (oe, stem, exts, exts.label="Time", fun=dsc, ret="oe", na.listwise=F, ...) {
  colNames <- ph.oe.ext(oe=oe, stem=stem, exts=exts)
  suboe <- oe[colNames]
  if (na.listwise==T)
    suboe <- na.omit(suboe)
  res <- sapply(suboe, fun, ...)
  if (all(sapply(res, is.null))) {
    invisible(NULL)
  } else {
    if (ret == "oe") {
      if (all(is.atomic(res)) & !is.null(names(res)) & length(res) == length(exts)) {
        retoe <- data.frame(exts, res, stringsAsFactors=F)
        names(retoe) <- c(exts.label, stem)
        rownames(retoe) <- names(res)
        return (retoe)
      }
    } else {
      res
    }
  }
}

# @returns any stems that have a certain suffix
ph.oe.ext.getStems <- function (oe, ext, numericOnly=F) {
  if (numericOnly==T) {
    nms <- names(oe)[sapply(oe, is.numeric)==T]
  } else {
    nms <- names(oe)
  }
  nms <- nms[ph.strs.detect.suffix.vec(nms, ext)]
  return (strs.replace.vec(nms, ext, ""))
}

# @returns colnames that have any of the passed exts as suffix
ph.oe.ext.getExtedColnames <- function (oe, exts) {
  nms <- names(oe)
  nmsLogi <- logical(length(nms)) # starts as FALSE
  for (ext in exts) {
    nmsLogi <- nmsLogi | ph.strs.detect.suffix.vec(nms, ext)
  }
  nms[nmsLogi==T]
}


# Used for split-half reliability analyses, this basically does a pairwise deletion (actually, an "NAing" not deletion) for 2 oes.
# So, if one oe is NA in cell [5, 4], this makes the other oe NA in cell [5. 4] too.
# Returns a list of 2 oes.
ph.oe.matchNAsOfTwooes <- function (oe.1, oe.2) {
  if ( (nrow(oe.1) != nrow(oe.2)) | (ncol(oe.1) != ncol(oe.2)) )
    stop("oe.matchNAsToOtheroe(): the 2 oes are different size, so I can't match NAs for them.")
  for (iColNum in 1:ncol(oe.1)) {
    oe.1[[iColNum]][is.na(oe.2[[iColNum]])] <- NA
    oe.2[[iColNum]][is.na(oe.1[[iColNum]])] <- NA
  }
  list(oe.1, oe.2)
}

# Moves column to the end of data frame (returns data frame with the column moved)
ph.oe.move.colToEnd <- function (oe, colName) {
  tmpCol <- oe[[colName]]
  oe[[colName]] <- NULL
  oe[[colName]] <- tmpCol
  oe
}

# Finds any columns that are all the same value and moves them to the last column positions.
ph.oe.move.uniformColsToEnd <- function (oe) {
  for (iColName in colnames(oe)) {
    if (all(is.na(oe[[iColName]]))==T)
      oe <- ph.oe.move.colToEnd(oe, iColName)
    if (!is.na(all(oe[[iColName]] == oe[[1, iColName]])) & all(oe[[iColName]] == oe[[1, iColName]])==T)
      oe <- ph.oe.move.colToEnd(oe, iColName)
  }
  oe
}

# Converts numeric cols to integer, if all values are equal to integers
# Helpful to run before outputting to SPSS
ph.oe.convert.numericToInteger <- function (oe) {
  for (iCol in 1:ncol(oe)) {
    if (class(oe[[iCol]])[1] == 'numeric') {
      col.noNAs <- na.omit(oe[[iCol]])
      if (all(col.noNAs == as.integer(col.noNAs)))
        oe[[iCol]] <- as.integer(oe[[iCol]])
    }
  }
  oe
}

# Prints a data frame proper, without column names above
# Note - this is a total hack. I replace names(oe) with the contents of row 1, completely invalidating the data structure,
#        which is why I just print it and do not return anything.
ph.oe.print <- function (oe, digits = NULL, right = TRUE, hideColNames = TRUE, row.names = FALSE, ...) {
  if (is.null(digits))
    digits <- getOption("digits")
  if (hideColNames==T) {
    lst <- as.list(oe[1, ])
    lst <- lapply(lst, function (item) { ifelse (is.numeric(item), format(item, digits=digits), as.character(item)) } )
    newNames <- sapply(lst, as.character)
    for (i in seq_along(newNames)) {
      if (nchar(newNames[i]) < nchar(names(oe)[i]))
        newNames[i] <- ph.str.padToLength(newNames[i], model=names(oe)[i])
    }
    names(oe) <- newNames
    oe <- oe[-1, ]
  }
  print(oe, digits=digits, right=right, row.names=row.names, ...)
  invisible(NULL)
}


# Removes columns from a data frame... you specify what to keep and/or remove
# Param keepOnly lets you just keep column names that *contain* any of the passed strings like keepOnly=c("threat", "p.")
# Param remove lets you remove any column names that *contain* any of the passed strings, even if they matched keepOnly
ph.oe.remove.cols <- function (oe, remove=NULL, remove.ifPrefix=NULL, remove.ifAllSame=F, remove.ifType=NULL, keepOnly=NULL, keepOnly.ifPrefix=NULL, keepOnly.ifType=NULL, silent=T) {
  nColsBefore <- ncol(oe)
  if (nColsBefore == 0)
    return (oe)
  colNumsToKeep <- 1:ncol(oe)
  if (is.null(keepOnly) & is.null(keepOnly.ifPrefix) & is.null(keepOnly.ifType)) # implies should keep all
    colNumsToKeep <- 1:ncol(oe)
  else
    colNumsToKeep <- NULL
  if (!is.null(keepOnly)) { 
    for (keepNeedle in keepOnly) {
      colNumsToKeep <- union(colNumsToKeep, which(ph.strs.detect.vec(names(oe), keepNeedle)))
    }
  }
  if (!is.null(keepOnly.ifPrefix)) {
    for (keepNeedle.prefix in keepOnly.ifPrefix) {
      colNumsToKeep <- union(colNumsToKeep, which(ph.strs.detect.prefix.vec(names(oe), keepNeedle.prefix)))
    }
  }
  if (!is.null(keepOnly.ifType)) {
    for (keepNeedle.type in keepOnly.ifType) {
      for (iColNum in 1:ncol(oe)) {
        if (is(oe[0, iColNum], keepNeedle.type)==T)
          colNumsToKeep <- union(colNumsToKeep, iColNum)
      }
    }
  }
  colNumsToRemove <- NULL
  if (!is.null(remove)) {
    for (removeNeedle in remove) {
      colNumsToRemove <- union(colNumsToRemove, which(ph.strs.detect.vec(names(oe), removeNeedle)))
    }
  }
  if (!is.null(remove.ifPrefix)) {
    for (removeNeedle in remove.ifPrefix) {
      colNumsToRemove <- union(colNumsToRemove, which(ph.strs.detect.prefix.vec(names(oe), removeNeedle)))
    }
  }
  if (!is.null(remove.ifType)) {
    for (removeNeedle.type in remove.ifType) {
      for (iColNum in 1:ncol(oe)) {
        if (is(oe[0, iColNum], removeNeedle.type)==T)
          colNumsToRemove <- union(colNumsToRemove, iColNum)
      }
    }
  }
  if (remove.ifAllSame==T) {
    for (iColNum in 1:ncol(oe)) {
      if (all(is.na(oe[[iColNum]]))==T)
        colNumsToRemove <- union(colNumsToRemove, iColNum)
      if (!is.na(all(oe[[iColNum]] == oe[[1, iColNum]])) & all(oe[[iColNum]] == oe[[1, iColNum]])==T)
        colNumsToRemove <- union(colNumsToRemove, iColNum)
    }
  }
  if (length(colNumsToKeep) > 0) {
    for (iColNum in colNumsToKeep) {
      if (iColNum %in% colNumsToRemove)
        colNumsToKeep <- colNumsToKeep[colNumsToKeep != iColNum]
    }
  }
  if (silent==F)
    ca("ph.oe.remove.cols(): Removed ", nColsBefore - length(unique(colNumsToKeep)), " cols, ", length(unique(colNumsToKeep)), " remaining\n")
  oe[colNumsToKeep]
}



# Recodes a vector, changing instances of the left sides of your pairs to whatever's on the right side
# Uses the recode() function from the package "car"... point is, this outputs a list describing the changes you're making, for your records--like to email to collaborators
# See car::recode() help for other syntaxes for the recodes param, other cool ways to specify your recodes
# @returns a list with:
#   $vec - the recoded vector
#   $records - a vector of strings which are the records of the recodes
#   $valueLabels - the new value labels, if any were passed in
# Usage e.g., change all occurrences in playerNames vec from regular names to nicknames:
# nickNamesL <- ph.vec.recode(vec = playerNames,
#                          recodes = "'Will' = 'Meteos';
#                                     'Zachary' = 'SneakyCastro'",
#                          valueLabels = c('Will' = 'Cloud9 team jungler',
#                                          'Zachary' = 'Cloud9 team fighter'))
# IMPORTANT NOTE - the syntax for recodes is different than for valueLabels. This is because I programmed the valueLabels one, whereas recodes comes from car.
ph.vec.recode <- function (vec, # vector of data to be recoded
                        recodes, # a named vector conveying your desired recodes... it's like assignment, left side becomes right
                        valueLabels=NULL, # named vector of labels describing the *old* values, such as from read.SPSS(). So in this one, each pair should be "value"="description of value"
                        silent=F
) {
  oldValues <- names(recodes)
  if (!all(oldValues %in% vec))
    stop("ph.vec.recode(): couldn't find all old values (left side of your recodes named vector) in vec")
  newVec <- car::recode(var=vec, recodes=recodes)
  u.old <- unique(vec)
  u.new <- unique(newVec)
  if (length(u.old) != length(u.new))
    opm("ph.vec.recode(): WARNING: the following rename records will display incorrectly, as you did a non 1-to-1 recode (the actual recode will still work")
  recordsVec <- NULL
  for (oldVal in u.old) {
    if (!(oldVal %in% u.new)) {
      newVal <- u.new[u.old == oldVal]
      recordStr <- "Renamed '"%_%oldVal%_%"' --> '" %_% newVal %_% "'"
      if (!is.null(valueLabels)) {
        recordStr <- recordStr %_% ", label: '" %_% valueLabels[[oldVal]] %_% "'"
        names(valueLabels)[names(valueLabels) == oldVal] <- newVal
      }
      if (silent==F)
        opm(recordStr)
      recordsVec <- c(recordsVec, recordStr)
    }
  }
  return (list(vec = newVec, records = recordsVec, valueLabels = valueLabels))
}



# @returns indexes of outliers
# @returns vec of indexes of outliers
# MAD is based on MAD-Median rule in Wilcox (2012, red book "Modern Statistics"), Section 3.13.4, p.97
#   For MAD, param is number of MADs an element must be to be deemed an outlier. So, to get like data-entry-error type outliers, use 15. For few outliers, 5.
ph.vec.outl.indexes <- function (vec, method="MAD", param=NULL) {
  if (!is.numeric(vec))
    return (integer()) # Say no outliers if non numeric
  switch(method,
         "MAD" = {
           if (is.null(param))
             param <- 2.24
           if (mad(vec, na.rm=T) == 0)
             return (integer())
           WRS::out(vec, crit=param)$out.id
         },
{
  stop("outliers.vec(): Invalid method: '"%_%method%_%"'")
})
}

# @returns logical vector, same length as input vec, with outlier positions as TRUE, non-outliers FALSE
ph.vec.outl.logi <- function (vec, req.cases=3, method="MAD", param=NULL) {
  if (sum(!is.na(vec)) < req.cases) {
    logical(length(vec))
  } else {
    indexes <- ph.vec.outl.indexes(vec=vec, method=method, param=param)
    seq_along(vec) %in% indexes
  }
}

# @returns data frame, same frame as input oe, and keeping cn.id same, but all other data are replaced by TRUE for outlier, FALSE for non-outlier
# @req.cases is minimum number of non-NA cases needed in the column to declare anything in that column an outlier
ph.oe.outl.logi <- function (oe, cn.id=names(oe)[1], cols=names(oe), req.cases=3, method="MAD", param=NULL) {
  origFirst <- oe[[cn.id]]
  outloe <- as.data.frame(lapply(oe, FUN=vec.outl.logi, req.cases=req.cases, method=method, param=param))
  outloe[[cn.id]] <- origFirst
  outloe
}

# @returns list of columns that are outliers for the passed row id.
ph.oe.outl.row  <- function (id, oe=data.frame(), cn.id=names(oe)[1], cols=names(oe), req.cases=3, method="MAD", param=NULL, outloe.premade=data.frame()) {
  if (identical(outloe.premade, data.frame())) {
    if (identical(oe, data.frame()))
      stop("oe.outl.row(): You didn't pass in a data frame.")
    outloe <- ph.oe.outl.logi(oe=oe, cn.id=cn.id, cols=cols, req.cases=req.cases, method=method, param=param)
  } else {
    outloe <- outloe.premade
  }
  if (sum(outloe[[cn.id]] == id) != 1)
    stop("oe.outl.row(): id must match exactly 1 row")
  rowVec <- as.vector(t(outloe[outloe[[cn.id]] == id, names(outloe) != cn.id]))
  names(outloe[names(outloe) != cn.id])[rowVec == T]
}

# @returns a 3-column oe, with the id column as col1, number of outliers in the row in col2, and list of which vars this row was outlier on in col3
# May want to make this just work on ONE row... returns a list, including a vector of colnames where this row has outliers
ph.oe.outl.rows <- function (oe, cn.id=names(oe)[1], cols=names(oe), req.cases=3, maxColNamesToShow=5, method="MAD", param=NULL, show.allRows=F) {
  outloe <- ph.oe.outl.logi(oe=oe, cn.id=cn.id, cols=cols, req.cases=req.cases, method=method, param=param)
  reportoe <- outloe[cn.id]
  reportoe$n <- apply(outloe[names(outloe) != cn.id], 1,
                      FUN=sum)
  
  priv.getColsString<-function (..., maxColNamesToShow) {
    colNames <- oe.outl.row(...)
    if (length(colNames) > maxColNamesToShow) {
      "'" %_% pa(colNames[1:maxColNamesToShow], collapse="', '") %_%
        "' ... ("%_%(length(colNames) - maxColNamesToShow)%_%" more)"
    } else if (length(colNames) > 0) {
      "'" %_% pa(colNames, collapse="', '") %_% "'"
    } else {
      character()
    }
  }
  
  reportoe$cols <- sapply(outloe[[cn.id]],
                          FUN=priv.getColsString,
                          cn.id=cn.id, cols=cols, method=method, param=param, maxColNamesToShow=maxColNamesToShow,
                          outloe.premade=outloe)
  if (show.allRows==T) {
    reportoe
  }
  else {
    reportoe[reportoe$n > 0, ]
  }
}



# Makes oe for output... very simple, just to make the whole outoe and rowoe rigmarole slightly easier
# (Not sure if this actually helps much)
ph.outoe.make <- function (colNames) {
  assign("OUTOE.SCRATCH", ph.data.frame.empty(colNames), envir=baseenv())
  assign("OUTOE.COLNAMES", colNames, envir=baseenv())
}
ph.outoe.push <- function (row) {
  oe <- get("OUTOE.SCRATCH", envir=baseenv())
  if (is.vector(row)) {
    rowoe <- data.frame(t(row))
    names(rowoe) <- names(oe)
    oe <- rbind(oe, rowoe)
  } else if (is.data.frame(row)) {
    names(row) <- names(oe)
    oe <- rbind(oe, row)
  }
  assign("OUTOE.SCRATCH", oe, envir=baseenv())
}
ph.outoe.get <- function () {
  outoe <- get("OUTOE.SCRATCH", envir=baseenv())
  names(outoe) <- get("OUTOE.COLNAMES", envir=baseenv())
  outoe
}

ph.str.padToLength <- function (stringToPad, model=NA, targetLength=NA, left=T, padChar=" ") {
  if (is.na(model) & is.na(targetLength))
    stop("ph.str.padToLength(): You didn't properly tell me what the targetLength you want is.")
  if (!is.na(model))
    targetLength <- nchar(model)
  if (targetLength <= nchar(stringToPad))
    return (stringToPad)
  pa(
    pa(rep(padChar, targetLength - nchar(stringToPad)), collapse=""),
    stringToPad)
}

# paste convenience aliases
pa <- function(..., sep="") paste(..., sep=sep)
ph.pa <- pa
"%_%" <- pa
pal <- function(...) pa(..., "\n")
ph.pal <- pal

# Takes a filename or path such as C:/Dropbox/aoe.csv and returns a named character vector with whole, path (just the directory part), core, and ext
parse.fname <- function (whole) {
  whole <- str_replace_all(whole, "/", "\\\\")
  if (ph.strs.detect(whole, ".")) {
    locMatrix <- str_locate_all(whole, "\\.")[[1]]
    lastDotLoc <- locMatrix[[nrow(locMatrix), 1]]
    ext <- substr(whole, lastDotLoc + 1, nchar(whole))
  } else
    ext <- ""
  if (ph.strs.detect(whole, "\\")) {
    locMatrix <- str_locate_all(whole, "\\\\")[[1]]
    lastBackslashLoc <- locMatrix[[nrow(locMatrix), 1]]
    path <- substr(whole, 1, lastBackslashLoc - 1)
  } else
    path <- "."
  core <- substr(whole,
                 ifelse(path == ".", 1, nchar(path) + 2),
                 ifelse(ext == "", nchar(whole), nchar(whole) - nchar(ext) - 1))
  fname.parsed <- parse.fname.makeFromParts(c(path=path, core=core, ext=ext))
  if (whole != fname.parsed[["whole"]])
    stop(pa("parse.fname(): parsing failed... the recombined whole, \"", fname.parsed[["whole"]], "\", does not match the original whole, \"", whole, "\""))
  fname.parsed
}
# Makes or remakes the "whole" element... you can pass it a complete fname.parsed or just the parts, i.e. a named char vector with path, core, and ext
parse.fname.makeFromParts <- function (fname.parts) {
  if (fname.parts[["path"]] == "")
    stop("parse.fname.remakeWhole(): fname.parsed[[\"path\"]] can't be empty... you should make it \".\"")
  fname.parts[["proper"]] <- pa(fname.parts[["core"]], ifelse(fname.parts[["ext"]] == "", "", pa(".", fname.parts[["ext"]])))
  fname.parts[["whole"]] <- pa(ifelse(fname.parts[["path"]] == ".", "", pa(fname.parts[["path"]], "\\")),
                               fname.parts[["core"]],
                               ifelse(fname.parts[["ext"]] == "", "", pa(".", fname.parts[["ext"]])))
  fname.parts
}

# Record processing steps to show number and percentage lost at each step
ph.recordStep <- function (nBefore, nAfter, desc, silent=F) {
  if (ph.isMode.skipOutput()) return (invisible(NULL))
  if (!exists("recordStepFile")) return ()
  nLost <- nBefore - nAfter
  pctLost <- 1 - (nAfter / nBefore)
  write(paste(sep="", nBefore, ",", 
              nAfter, ",",
              nLost, ",",
              pctLost * 100, ",",
              desc, ","),
        file=recordStepFile, append=TRUE)
  if (silent==F)
    ca(sep="", "Recording step: ", nBefore, " -> ", nAfter, ", ", nLost, " (", format(pctLost * 100, digits=3), "%) lost. Desc: \"", desc, "\"\n")
}

# Reads entire file into character string, like PHP's file_get_contents() function.
ph.read.entireFile <- function (fname) {
  fconn <- file(fname, "r")
  len <- file.info(fname)$size
  toret <- readChar(fconn, len)
  close(fconn)
  toret
}

# Init record processing steps
# Pass filename, or pass "" to stop recording
ph.recordStepInit <- function (filename) {
  if (ph.isMode.skipOutput()) return (invisible(NULL))
  if (filename == "") {
    if (exists("recordStepFile"))
      rm(recordStepFile, envir=globalenv())
    return ()
  }
  recordStepFile <<- filename
  write("nBefore,nAfter,nLost,pctLost,desc", filename)
}

# Returns a random alphanumeric string of specified length
ph.rndString <- function (len) {
  paste(sample(c(letters, 0:9), len, replace=T), collapse="")
}

# Does a shell, replacing slashes with double-backslashes, to make it work
ph.run.shell <- function (cmd, ...) {
  shell(str_replace_all(cmd, "/", "\\\\"), ...)
}

# Does a shell.exec, replacing slashes with double-backslashes, because if you don't do that then relative paths that start with a folder name, shell.exec("opn/file.txt"), don't work 
ph.run.file <- function (path, ...) {
  if (ph.isMode.skipOutput()) return (invisible(NULL))
  shell.exec(str_replace_all(path, "/", "\\\\"), ...)
}

# Creates scratch variables for sandbox play
ph.sandbox <- function () {
  salist <<- list(5, 10, 12)
  soe <<- data.frame(i=c(3, 5, 7), j=c(9, 10, 200), f=c("item1", NA, "item3"), l=c(TRUE, TRUE, FALSE), s=c("hi", "there", "how are you"))
  cal("SANDBOX... Created: ")
  cat ("salist: \n"); print(salist); cal ("data frame soe: ");print(soe)
}

# Calls str() with just the top 2 levels
ph.strtop <- function (...) {
  str(..., max.level=2)
}


# Returns start position of (first occurrence of) string needle in haystack, NOT using regex
# Returns FALSE if not found
ph.strs.pos <- function (haystack, needle, caseSens=T, findLastOne=F, na.returnsFalse=T) {
  if (is.null(haystack)) {
    stop("ph.strs.pos(): Your haystack arg is null -Phil")
  } else if (is.na(haystack)) {
    if (na.returnsFalse==T) {
      return (FALSE)
    } else {
      stop("ph.strs.pos(): Your haystack arg is NA -Phil")
    }
  } else if (length(haystack) != 1) {
    stop("ph.strs.pos(): Your haystack arg is length != 1 and should be just 1 -Phil")
  } else if (is.null(needle)) {
    stop("ph.strs.pos(): Your needle arg is null -Phil")
  } else if (is.na(needle)) {
    if (na.returnsFalse==T) {
      return (FALSE)
    } else {
      stop("ph.strs.pos(): Your needle arg is NA -Phil")
    }
  } else if (length(needle) != 1) {
    stop("ph.strs.pos(): Your needle arg is length != 1 and should be just 1 -Phil")
  } else if (nchar(haystack) < nchar(needle)) {
    return (FALSE)
  } else if (needle == "") {
    return (ifelse(haystack == "", 1, FALSE))
  }
  if (caseSens==F) {
    haystack <- tolower(haystack)
    needle <- tolower(needle)
  }
  res <- gregexpr(pattern=needle, text=haystack, fixed=T)
  if (res[[1]][1] == -1) {
    FALSE
  } else {
    if (findLastOne==F) {
      res[[1]][1]
    } else {
      res[[1]][length(res[[1]])]
    }
  }
}
ph.strs.pos.vec <- function (vec, needle, caseSens=T) {
  if (length(vec) < 1)
    stop("Your vec arg is 0-length -Phil")
  posVec <- vector()
  sapply( vec,
          FUN=function (applyItem)
          { ph.strs.pos(applyItem, needle, caseSens) } )
}

ph.strs.detect <- function (haystack, needle, caseSens=T) {
  ph.strs.pos(haystack, needle, caseSens=T) != FALSE
}
ph.strs.detect.vec <- function (vec, needle, caseSens=T) {
  ph.strs.pos.vec(vec, needle, caseSens) != 0
}

ph.strs.detect.prefix <- function (haystack, needle, caseSens=T) {
  ph.strs.pos(substr(haystack, 1, nchar(needle)), needle, caseSens) != FALSE
}
ph.strs.detect.prefix.vec <- function (vec, needle, caseSens=T) {
  ph.strs.detect.vec(substr(vec, 1, nchar(needle)), needle, caseSens)
}

ph.strs.detect.suffix <- function (haystack, needle, caseSens=T) {
  ph.strs.pos(substr(haystack, nchar(haystack)+1 - nchar(needle), nchar(haystack)), needle, caseSens) != FALSE
}
ph.strs.detect.suffix.vec <- function (vec, needle, caseSens=T) {
  if (length(vec) < 1)
    stop("Your vec arg is 0-length -Phil")
  sapply( vec,
          FUN=function (applyItem)
          { ph.strs.detect.suffix(applyItem, needle, caseSens) } )
}

# Returns the part of haystack that's before needle
# By default, it takes the first hit for needle. Or, set fromLast=T, and it'll extract after the last hit.
ph.strs.extractBefore <- function (haystack, needle, fromLast=F) {
  substr(haystack, 1,
         ph.strs.pos(haystack, needle, findLastOne=fromLast) - 1)
}
# Vector version of ph.strs.extractBefore()
ph.strs.extractBefore.vec <- function (vec, needle, fromLast=F) {
  if (length(vec) < 1)
    stop("Your vec arg is zero-length -Phil")
  sapply( vec,
          FUN=function (applyItem)
          { ph.strs.extractBefore(applyItem, needle, fromLast=fromLast) } )
}

# Returns the substring of passed string between the two passed strings (usually chars) non-inclusive
#   Returns "" if there's nothing between or FALSE if needle1 or needle2 isn't found
#   Setting inclusive=T will mean the return includes needle1 and needle2 on the ends
ph.strs.extractBetween <- function (haystack, needle1, needle2, inclusive=F) {
  needle1Pos <- ph.strs.pos(haystack, needle1)
  needle2Pos <- ph.strs.pos(haystack, needle2)
  if (needle1Pos == F | needle2Pos == F) {
    FALSE
  } else {
    haystack.afterNeedle1 <- substr(haystack, needle1Pos + nchar(needle1), nchar(haystack))
    endPos.within.haystack.afterNeedle1 <- ph.strs.pos(haystack.afterNeedle1, needle2)
    if (endPos.within.haystack.afterNeedle1 == F) {
      return (FALSE)
    } else {
      if (inclusive==F) {
        substr(haystack.afterNeedle1, 1, endPos.within.haystack.afterNeedle1-1)
      } else {
        pa(needle1, substr(haystack.afterNeedle1, 1, endPos.within.haystack.afterNeedle1-1), needle2)
      }
      
    }
  }
}
# Vector version of ph.strs.extractBetween()
ph.strs.extractBetween.vec <- function (vec, needle1, needle2, inclusive=F) {
  if (length(vec) < 1)
    stop("Your vec arg is zero-length -Phil")
  sapply( vec,
          FUN=function (applyItem)
          { ph.strs.extractBetween(applyItem, needle1, needle2, inclusive) } )
}

# Returns the rest of haystack after needle
# By default, it takes the first hit for needle. Or, set fromLast=T, and it'll extract after the last hit.
ph.strs.extractAfter <- function (haystack, needle, fromLast=F) {
  substr(haystack,
         ph.strs.pos(haystack, needle, findLastOne=fromLast) + nchar(needle),
         nchar(haystack))
}
# Vector version of ph.strs.extractAfter()
ph.strs.extractAfter.vec <- function (vec, needle, fromLast=F) {
  if (length(vec) < 1)
    stop("Your vec arg is zero-length -Phil")
  sapply( vec,
          FUN=function (applyItem)
          { ph.strs.extractAfter(applyItem, needle, fromLast=fromLast) } )
}

# Removes the characters between the passed needle chars/strings.
# Returns whole haystack if needles aren't found.
ph.strs.removeBetween <- function (haystack, needle1, needle2, inclusive=F) {
  ph.strs.replace(haystack,
               ph.strs.extractBetween(haystack, needle1, needle2, inclusive=inclusive),
               "")
}
# Vector version of ph.strs.removeBetween()
ph.strs.removeBetween.vec <- function (vec, needle1, needle2, inclusive=F) {
  if (length(vec) < 1)
    stop("Your vec arg is zero-length -Phil")
  sapply( vec,
          FUN=function (applyItem)
          { ph.strs.removeBetween(applyItem, needle1, needle2, inclusive) } )
}

ph.strs.replace <- function(haystack, needle.old, needle.new) {
  gsub(pattern=needle.old, replacement=needle.new, x=haystack, fixed=T)
}
ph.strs.replace.vec <- Vectorize(ph.strs.replace, USE.NAMES=F)

ph.strs.replace.suffix <- function(haystack, needle.old, needle.new, caseSens=T) {
  if(length(haystack) == 0) {
    haystack
  } else if (ph.strs.pos(substr(haystack, nchar(haystack)+1 - nchar(needle.old), nchar(haystack)), needle.old, caseSens) != FALSE) {
    return (pa(substr(haystack, 1, nchar(haystack) - nchar(needle.old)), needle.new))
  } else {
    haystack
  }
}
ph.strs.replace.suffix.vec <- Vectorize(ph.strs.replace.suffix, USE.NAMES=F)


# Returns passed string with all spaces, tabs, and newlines removed.
ph.str.stripSpaces <- function (text) {
  gsub("\\s","", text)
}

# Returns passed string with spaces (or passed strToRemove) removed from start or end of string only.
ph.str.stripAtStart <- function (text, strToRemove=" ") {
  hitMatrix <- str_locate_all(text, strToRemove)[[1]]
  if (nrow(hitMatrix) == 0)
    return (text)
  toret <- text; lastEnd <- 0 # loop prep
  for (iRow in 1:nrow(hitMatrix)) {
    if (hitMatrix[iRow, "start"] == lastEnd + 1) {
      toret <- str_replace(toret, strToRemove, "")
      lastEnd <- hitMatrix[iRow, "end"]
    }
  }
  toret
}

# Returns legal varname.
# You can use this to check if text is legal by checking nchar(ph.str.make.legal.varname(text)) == nchar(text)
ph.str.make.legal.varname <- function (name) {
  str_replace_all(name, "[- \\[\\$@\\\\/:\\*<>\\|\\\"]", "")
}
ph.str.is.legal.varname <- function (name) { ph.str.make.legal.varname(name) == name }

# Returns legal filename.
ph.str.make.legal.filename <- function (name) {
  name <- str_trim(name)
  name <- str_replace_all(name, "[\\[\\$@\\\\/:\\*<>\\|\\\"]", "'")
  name <- ph.str.stripAtStart(name, "'")
  str_replace_all(name, "[\\[\\$@\\\\/:\\*<>\\|\\\"]", "'")
}

# Removes the characters at the beginning of passed string if they match needle.
ph.str.removePrefix <- function (haystack, needle) {
  if (substr(haystack, 1, nchar(needle)) == needle)
    substr(haystack, nchar(needle)+1, nchar(haystack))
  else
    haystack
}

# Removes the characters at the beginning of passed string if they match needle.
ph.str.removePrefix.vec <- function (vec, needle) {
  if (length(vec) < 1)
    stop("Your vec arg is not legit -Phil")
  retVec <- vector()
  for (i in 1:length(vec))
    retVec[i] <- ph.str.removePrefix(vec[i], needle)
  retVec
}

# Passed a list, returns a list with n elements, each 
ph.splitListIntoEqualChunks <- function(dataList, n) {
  dataList <- as.list(dataList)
  sliceLength <- length(dataList) / n
  i <- seq_along(dataList)
  split(dataList, ceiling(i/sliceLength))
}

# Passed an oe, returns a list of oes, where each oe chunk has (approximately) 1/n rows
ph.splitoeIntoEqualChunks <- function(oe, n) {
  sliceLength <- nrow(oe) / n
  i <- 1:nrow(oe)
  split(oe, ceiling(i/sliceLength))
}

# Passed an oe, returns a list of oes, where each oe chunk has (approximately) 1/n rows
# - All rows of any given delimiter value will be in the same chunk, no breakups.
# - Each chunk is greedy, as in it can extend longer than sliceLength, because it looks below but not above aprioriCutoff to extend the chunk based on delim.
# - So keep in mind, this function may return fewer chunks than the passed n.
ph.splitoeIntoEqualChunksWithoutBreakingDelimiter <- function(oe, n, delim) {
  if (nrow(oe) == 0)
    stop(pa("ph.splitoeIntoEqualChunksWithoutBreakingDelimiter(): You asked me to split an empty oe into multiple chunks. -Phil"))
  if (nrow(oe) == 1)
    return (list(oe))
  sliceLength <- ceiling(nrow(oe) / n)
  chunkList <- list()
  chunkStart <- 1
  # ca("a priori sliceLength = ", sliceLength, "; Splitting into list without breaking delimiter: ", n, " chunks","\n")
  for (iChunk in 1:(n)) {
    aprioriCutoff <- iChunk * sliceLength
    if (aprioriCutoff <= chunkStart) next
    if (aprioriCutoff < nrow(oe)) {
      iDelim <- aprioriCutoff + 1 # hackish, to make aprioriCutoff the last element of the chunk
      if (delim[aprioriCutoff] == delim[aprioriCutoff+1]) {
        for (iDelim in (aprioriCutoff+1):nrow(oe)) {
          if (delim[iDelim] != delim[iDelim-1]) {
            break
          }
        }
      }
    } else {
      chunkList <- c(chunkList, list(oe[(chunkStart:nrow(oe)), ]))
      return (chunkList)
    }
    if (iDelim == nrow(oe)) { # means didn't find a change in delim, having searched the entire rest of delim
      chunkList <- c(chunkList, list(oe[(chunkStart:nrow(oe)), ]))
      return (chunkList)
    }
    if (chunkStart >= iDelim-1)
      stop(pa("ph.splitoeIntoEqualChunksWithoutBreakingDelimiter(): You asked me to split just \"", nrow(oe), "\" row into multiple chunks. -Phil"))
    chunkList <- c(chunkList, list(oe[chunkStart:(iDelim-1), ]))
    if (iDelim == nrow(oe)) {
      chunkList <- c(chunkList, list(oe[iDelim, ]))
      return (chunkList)
    }
    chunkStart <- iDelim
  }
  chunkList
}

tb <- traceback

ph.tim <- function (...) system.time(...)[3]

# Shows the timer
ph.timeCheck <- function (desc="") {
  if (desc != "") desc  <- pa(" : ", desc)
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  elapsedSec <- toc[[1]] - tic[[1]]
  if (elapsedSec >= 300) {
    pa(signif(elapsedSec / 60, digits=3), "min", desc)
  } else {
    pa(signif(elapsedSec, digits=3), "s", desc)
  }
}

# Starts a timer... flank some lines of code with timeIt and timeCheck so you can see how long it takes
ph.timeIt <- function(gcFirst=TRUE, type=c("elapsed", "user.self", "sys.self")) {
  type <- match.arg(type)
  assign(".type", type, envir=baseenv()) # he used the R base package namespace  to hold vars between functions, and so it doesn't interfere with my user namespace
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]         
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

# Starts USER's timer... flank some lines of code with timer.start and timer.check so you can see how long it takes
ph.timer.start <- function(gcFirst=TRUE, type=c("elapsed", "user.self", "sys.self")) {
  type <- match.arg(type)
  assign("usertimer.type", type, envir=baseenv()) # he used the R base package namespace  to hold vars between functions, and so it doesn't interfere with my user namespace
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]         
  assign("usertimer.tic", tic, envir=baseenv())
  invisible(tic)
}

# Shows the USER's timer
ph.timer.check <- function (desc="") {
  if (desc != "") desc  <- pa(" : ", desc)
  type <- get("usertimer.type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get("usertimer.tic", envir=baseenv())
  elapsedSec <- toc[[1]] - tic[[1]]
  if (elapsedSec >= 300) {
    pa(signif(elapsedSec / 60, digits=3), "min", desc)
  } else {
    pa(signif(elapsedSec, digits=3), "s", desc)
  }
}

# returns time stamp including time since last ph.timeIt()
ph.timeStamp <- function(time=T, timer=F) {
  if (time==T & timer==T)
    format(Sys.time(), pa("(%H:%M:%S|", ph.timeCheck(), ")"))
  else if (time==F & timer==F)
    ""
  else if (time==T)
    format(Sys.time(), pa("(%H:%M:%S|)"))
  else if (timer==T)
    pa("(", ph.timeCheck(), ")")
}

# returns string w/o leading or trailing whitespace.
# src: http://stackoverflow.com/questions/2261079/whitespace-in-r
ph.trim <- function (x) {
  gsub("^\\s+|\\s+$", "", x)
}

ph.vars.list <- function (varNames, nameColName="varName", valueColName="value") {
  oe <- data.frame()
  for (iVarName in varNames) {
    oe <- rbind(oe, data.frame(iVarName, get(iVarName)))
  }
  names(oe) <- c(nameColName, valueColName)
  oe
}

# @Returns TRUE if this vector is good for statistical tests... has to be numeric and have variance
ph.vec.is.goodForTests <- function (x) {
  if (is.null(x))
    return (FALSE)
  if (all(is.na(x)) == T)
    return (FALSE)
  if (!is.numeric(x))
    return (FALSE)
  if (sd(x, na.rm=T) == 0)
    return (FALSE)
  TRUE
}

# Removes elements from a vec
ph.vec.remove.elements <- function (vec, remove=NULL, remove.ifPrefix=NULL, keepOnly=NULL, keepOnly.ifPrefix=NULL, silent=T) {
  nElemsBefore <- length(vec)
  if (nElemsBefore == 0)
    vec
  indexesToKeep <- 1:length(vec)
  if (is.null(keepOnly) & is.null(keepOnly.ifPrefix)) # implies should keep all
    indexesToKeep <- 1:length(vec)
  else
    indexesToKeep <- NULL
  if (!is.null(keepOnly)) { 
    for (keepNeedle in keepOnly) {
      indexesToKeep <- union(indexesToKeep, which(ph.strs.detect.vec(vec, keepNeedle)))
    }
  }
  if (!is.null(keepOnly.ifPrefix)) {
    for (keepNeedle.prefix in keepOnly.ifPrefix) {
      indexesToKeep <- union(indexesToKeep, which(ph.strs.detect.prefix.vec(vec, keepNeedle.prefix)))
    }
  }
  indexesToRemove <- NULL
  if (!is.null(remove)) {
    for (removeNeedle in remove) {
      indexesToRemove <- union(indexesToRemove, which(ph.strs.detect.vec(vec, removeNeedle)))
    }
  }
  if (!is.null(remove.ifPrefix)) {
    for (removeNeedle in remove.ifPrefix) {
      indexesToRemove <- union(indexesToRemove, which(ph.strs.detect.prefix.vec(vec, removeNeedle)))
    }
  }
  if (length(indexesToKeep) > 0) {
    for (iIndex in indexesToKeep) {
      if (iIndex %in% indexesToRemove)
        indexesToKeep <- indexesToKeep[indexesToKeep != iIndex]
    }
  }
  if (silent==F)
    ca("ph.vec.remove.elements(): Removed ", nElemsBefore - length(unique(indexesToKeep)), " elements, ", length(unique(indexesToKeep)), " remaining\n")
  if (length(indexesToKeep) > 0)
    vec[indexesToKeep]
  else
    vec[NULL]
}

# A shortcut function: Returns same vector but without any NAs
ph.wo.na <- function (vec) {
  vec[!is.na(vec)]
}



# ------------------------------------ < Startup code executed upon sourcing this file > ------------------------------------
ph.timeIt()
if (PH$REMOVE.PH.PREFIX)
  ph.removePHPrefix()

# This is EnviHack, see note above on EnviHack
PH$LS.AFTER <- ls()
PH$LIST <- setdiff(PH$LS.AFTER, PH$LS.BEFORE)
PH.ENV <- new.env()
for (obj in PH$LIST) {
  assign(obj, get(obj), envir=PH.ENV)
}

PH$LS.AFTER <- PH$LS.BEFORE <- NULL
rm(list=PH$LIST)
rm(obj)
if ('PH.ENV' %in% search())
  detach(PH.ENV)
attach(PH.ENV)
