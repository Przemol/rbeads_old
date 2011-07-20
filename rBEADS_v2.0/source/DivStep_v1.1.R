# TODO: Add comment
# 
# Author: przemol
###############################################################################

DivStep <- function(sample, control, enriched_regions=NULL, genome=Celegans) {
	
		catTime("Divideing", e={	
			divstep <- sapply(names(sample), function(x) sample[[x]] / control[[x]])
			divstep <- RleList(divstep)
		})

		catTime("Scaling", e={
			aa <- median(median(divstep, na.rm=T))
			divstep <- divstep / aa
		})
	cat('INFO: Scaling coefficient = ', aa, '\n', sep='')
	
	return(round(divstep, 3))
}

catTime <- function(..., e=NULL, file="", gc=FALSE) {
	cat(..., "...", sep="", file=file, append=TRUE)
	cat("\t<", system.time(e, gcFirst=gc)[3], "s>\n", sep="", file=file, append=TRUE)	
}

