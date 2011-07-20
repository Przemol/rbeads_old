# TODO: Add comment
# 
# Author: przemol
###############################################################################


MappabilityCorrection <- function(GCnormTrack, mappabilityTrack=NULL, genome=Celegans) {
	
	
	catTime("Mappability correction", e={			
		# 0) Old way (multyply than divide by the same value in DIV step)
		#cov_GCmap <- cov * cmap.fd
		# 1) masking regions of mappability lower than 100 (mappability score=0.25) from precalculated mappability track
		#GCnormTrack[cmap.fd < 0.25] <- NA
		if (genome@provider_version == 'ce6') {	
			#names(GCnormTrack) <- seqnames(Celegans)[c(1:4,7,5,6)]
			#GCnormTrack[mappabilityTrack[c(1:4,7,5,6)] < 100] <- NA
					for(i in names(GCnormTrack)) {
						cat(.)
						GCnormTrack[[i]][mappabilityTrack[[i]] < 100] <- NA
					}
		} else if (genome@provider_version == 'hg19') {	
			for(i in names(GCnormTrack)) {
				cat(.)
				GCnormTrack[[i]][mappabilityTrack[[i]] < 100] <- NA
			}
		} else if (genome@provider_version == 'mm9') {
			for(i in names(GCnormTrack)) {
				cat(.)
				GCnormTrack[[i]][mappabilityTrack[[i]] < 100] <- NA
			}
		} else {
			warning('Unnkonown genome!')
		}
		# 2) masking regions within lowest 5% of all read values
	# cov_GCmap[cov.r <= quantile(as.vector(cov.r), 0.05)] <- NA			
	})

	return(GCnormTrack)
	
}

catTime <- function(..., e=NULL, file="", gc=FALSE) {
	cat(..., "...", sep="", file=file, append=TRUE)
	cat("\t<", system.time(e, gcFirst=gc)[3], "s>\n", sep="", file=file, append=TRUE)	
}