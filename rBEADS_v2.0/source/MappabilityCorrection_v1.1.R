# TODO: Add comment
# 
# Author: przemol
###############################################################################


MappabilityCorrection <- function(GCnormTrack, mappabilityTrack=NULL) {
	
	if (is.null(mappabilityTrack)) {
		catTime("Loading default mappability track", e={			
			varname <- load("/Users/przemol/Documents/workspace/RBeads/rBeads/precalculated/MappabilityCe6_v2.Rdata")
			assign("mappabilityTrack", get(varname))
		})
		
	}
	
	catTime("Mappability correction", e={			
		# 0) Old way (multyply than divide by the same value in DIV step)
	#cov_GCmap <- cov * cmap.fd
		# 1) masking regions of mappability lower than 100 (mappability score=0.25) from precalculated mappability track
		#GCnormTrack[cmap.fd < 0.25] <- NA
	names(GCnormTrack) <- seqnames(Celegans)[c(1:4,7,5,6)]
	GCnormTrack[mappabilityTrack[c(1:4,7,5,6)] < 100] <- NA
		# 2) masking regions within lowest 5% of all read values
	# cov_GCmap[cov.r <= quantile(as.vector(cov.r), 0.05)] <- NA			
	})

	return(GCnormTrack)
}

catTime <- function(..., e=NULL, file="", gc=FALSE) {
	cat(..., "...", sep="", file=file, append=TRUE)
	cat("\t<", system.time(e, gcFirst=gc)[3], "s>\n", sep="", file=file, append=TRUE)	
}