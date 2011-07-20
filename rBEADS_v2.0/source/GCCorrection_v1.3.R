# TODO: Add comment

# Author: przemol
###############################################################################
require(GenomicRanges)
require(rtracklayer)

GCCorrection <- function(ranges.raw, enriched_regions, nonMappableFilter, genome=Celegans, desc='', smoothing_spline=FALSE, cutoff=c(35, 160)) {
	
	#Mask out reads in enriched regions
	if (!is.null(enriched_regions)) {
		catTime("Mask out reads in enriched regions", e={
			ERoverlaps <- findOverlaps(ranges.raw, enriched_regions, select="first") 	
		})
		#INFO: percentage of reads in enriched regions
		cat("\tINFO: percentage of reads in non-enriched regions: ", round(100 * sum(is.na(ERoverlaps))/length(ranges.raw), 2), "%\n", sep='')
	}
	
	#Calculate input GC content from 200bp extended reads (ranges.raw) ##TOO LONG! 221.484s
	catTime("Calculate input GCcontent", e={
		GCcontent <- as.integer(letterFrequency(getSeq(genome, ranges.raw, as.character=FALSE), "GC"))
	})
	
	#Sample genome for GCcontent
	catTime("Sample genome for GCcontent", e={
		if (!is.null(enriched_regions)) {
			#Calculate ogical vectors of non-enriched regions
			nonEnrichedRegionsLogi <- !coverage(GRanges(space(enriched_regions), unlist(ranges(enriched_regions)), "*", seqlengths=seqlengths(genome)[seqlevels(ranges.raw)]))
			#Perform logical sum of non-enriched regions and mappable regions
			nonEnrichedMappableRegionsLogi <- nonEnrichedRegionsLogi & nonMappableFilter[names(nonEnrichedRegionsLogi)]
		} else {
			nonEnrichedMappableRegionsLogi <- nonMappableFilter
		}
		#Calculate GC pecrentage among the chromosomes
		GC <- lapply(seqlevels(ranges.raw), function(x) hist( c(letterFrequencyInSlidingView(getSeq(genome, x, as.character=F), 200, "GC"), rep(NA, 199))[as.logical(nonEnrichedMappableRegionsLogi[[x]])], 0:200, plot=F ) ) #Bottle neck, long, memory consuming
		a.dens <- apply( sapply(GC, function(x) x$density), 1, function(x) weighted.mean(x, sum(nonEnrichedMappableRegionsLogi)[seqlevels(ranges.raw)]))
	})
	
	#INFO: percentage of genome to be sampled
	cat("\tINFO: percentage of genome to be sampled: ", sum(as.numeric(sum(nonEnrichedMappableRegionsLogi))) / sum(as.numeric(seqlengths(genome))), "\n") 

	#Calculate histograms for genomic (a) nad and sample (b) GC content
	catTime("Calculate histograms for genomic (a) nad and sample (b) GC content", e={
		if (!is.null(enriched_regions)) {
			b <- hist(GCcontent[is.na(ERoverlaps)], 0:200, plot=F)
		} else {
			b <- hist(GCcontent, 0:200, plot=F)
		}
	
		pdf(sprintf("IMG - GCdistribution - %s.pdf", desc), width = 12.0, height = 7.5, onefile = FALSE, paper = "special", encoding = "TeXtext.enc")
		plot(b$mids, b$density, col="red", type="l", main=sprintf("GCdistribution - %s", desc))
		lines(b$mids, a.dens, type="l", col="blue")
		legend("topleft", c("Non enriched GENOMIC GC content distribution", "Non enriched SAMPLE GC content distribution"), fill=c("blue", "red") )
		dev.off()
	})
	
	cat('\tINFO: scales <- c(', paste(a.dens/b$density, collapse=','), ')\n', sep='')
	
	#Calculate GC weighting vector
	catTime("Calculate GC weighting vector", e={												
		scales <- a.dens/b$density
	
		pdf(file=sprintf("IMG - GCoutlier - %s.pdf", desc), width = 12.0, height = 7.5, onefile = FALSE, paper = "special", encoding = "TeXtext.enc")
		plot(scales, main=sprintf("GCoutlier - %s", desc))
		
		#Remove outliers from the vecor
		#Grubbs test
		#while (grubbs.test(scales[,2])$p < 0.1 ) {
		#	scales <- scales[-which.max(scales[,2]), ]
		#}
		#Alternative ways of doing that
		#1) Quantile
		#scales <- scales[which(scales[,2] < quantile(scales[,2], .9)), ]
		#2) Boxplot
		scales[ scales %in% boxplot(scales, plot=F)$out ] <- NA
		#3) Fixed
		scales[c( 1:(cutoff[1]-1), (cutoff[2]+1):200 )] <- NA
		abline(h=0, v=c(cutoff[1], cutoff[2]), col = "gray60")
		
		
		#Roboust prediction of scaling ceficient by fitting a cubic smoothing spline to the supplied data.
		SSpoints <- cbind(1:200, scales)[!is.na(scales) & !(scales %in% boxplot(scales, plot=F)$out), ]	
		SS <- smooth.spline(SSpoints, df=10)
		scales.fit <- cbind(predict(SS, 0:200)$x, predict(SS, 0:200)$y)
		lines(scales.fit, col="green")
		
		if (smoothing_spline) {
			scales.f <- scales.fit
		}
		points(scales, col="red", pch="*")
		points(SSpoints, col="green", pch="o")
		dev.off()
	})
	
	#Assign GC score to every read
	catTime("Assign GC score to every read", e={
		GCcontent[GCcontent == 0] <- 1
		scores <- scales[GCcontent]
	})
	
	#Prepare GCscore-enriched GRanges
	catTime("Prepare GCscore-enriched GRanges", e={				
		values(ranges.raw) <- scores
		ranges.f1 <- ranges.raw[!is.na(values(ranges.raw)),]
	})
	
	#calculate GC weighted coverage with given accuracy
	catTime("calculate GC weighted coverage with given accuracy", e={
		acc = 1000
		w <- sapply( seqlevels(ranges.f1), function(x) values( ranges.f1[seqnames(ranges.f1) == x] )$value * acc )
		#Output coverage
		cov <- coverage(ranges.f1, weight=w) / acc
		cov.r <- round(cov)
	})

	#Mask nonGC correctable regions TOO DO: too log, too much memory
	catTime("Masking non-GCcorrectable regions", e={
		for(i in names(cov.r)) {
			cat('.')
			GCchr <- letterFrequencyInSlidingView(getSeq(genome, names=i, as.character=F, strand="*"), 200, "GC")
			notGCcorrectableRegions <- GRanges(seqnames=i, ranges=IRanges(which( ! (GCchr >= cutoff[1] & GCchr <= cutoff[2]) ), width=200))
			seqlengths(notGCcorrectableRegions) <- seqlengths(genome)[ i ]
			notGCcorrectableRegions <- c(notGCcorrectableRegions, GRanges(seqnames=i, ranges=IRanges( seqlengths(genome)[i]-199, width=200)) )
			cov.r[[i]][ (coverage(notGCcorrectableRegions) > 0)[[i]] ] <- NA
		}
	})

	return(cov.r)
}

catTime <- function(..., e=NULL, file="", gc=FALSE) {
	cat(..., "...", sep="", file=file, append=TRUE)
	cat("\t<", system.time(e, gcFirst=gc)[3], "s>\n", sep="", file=file, append=TRUE)	
}
