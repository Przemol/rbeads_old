# TODO: Add comment
# 
# Author: przemol
###############################################################################


# TODO: Add comment
# 
# Author: przemol
###############################################################################


require(GenomicRanges)
require(rtracklayer)
library(multicore)

EnrichedRegions <- function(ranges.raw, desc="EnrichedRegions1", genome=Celegans ) {
	
	
	cat("GC correction - peak calling", "\n")
	
	#INTEGRATED PEAK CALLER
	
	#Calculate coverage
	catTime("Calculate coverage", e={
				combExtCoverRep1 <- coverage(ranges.raw)
			})	
	
	#Do peak calling
	catTime("Do peak calling", e={
				a = quantile(quantile(combExtCoverRep1, 0.75), 0.75)
				enriched_regions <- slice(combExtCoverRep1, lower = a)
				peakSumsRep1 <-viewSums(enriched_regions)
				enriched_regions <- enriched_regions[peakSumsRep1 >= quantile(peakSumsRep1[[3]], .90)]
			})
	cat('INFO: ', 100 *sum(sum(width(enriched_regions))) / sum(as.numeric(seqlengths(genome)[seqlevels(ranges.raw)])), '% of genome in ER\n', sep='')
	#INFO: prepare the peak calling bed file to view in IGB/IGV
	catTime("INFO: prepare the peak calling bed file to view in IGB/IGV", e={															
				export.bed(enriched_regions, sprintf("%s.ERpeakCall.bed", desc))	
			})
	return(RangedData(enriched_regions))
}


