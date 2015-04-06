



# Functions, just names
ph.an.oe.3D.ttests <- ph.oea.3D.ttests
ph.an.oe.3D.to2D <-ph.oea.3D.to2D
ph.an.cor.set <- ph.oea.cor.set
ph.an.t.1sample.set <- ph.oea.t.1sample.set
ph.an.t.set <- ph.oea.t.set



# Functions, with different specifications
ph.an.coh.trt.s <- function (data, varList, preSuffix, postSuffix, grp1, grp2) {
  ph.oea.coh.trt.set(data=data, stems=varList, preSuffix=preSuffix, postSuffix=postSuffix) # Legacy
}

