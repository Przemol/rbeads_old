
# TODO: Add comment
# 
#frame_files <- lapply(sys.frames(), function(x) x$ofile)
#frame_files <- Filter(Negate(is.null), frame_files)
#PATH <- dirname(frame_files[[length(frame_files)]])

# Author: przemol
###############################################################################

require(Rsamtools)
require(BSgenome.Celegans.UCSC.ce6)
require(BSgenome.Hsapiens.UCSC.hg19)
require(BSgenome.Mmusculus.UCSC.mm9)

ImportBAM <- function(bam.file=dir(pattern="\\.bam$")[1], genome=Celegans, desc=unlist(strsplit(bam.file, "\\."))[1], resize_length=200, quality_cutoff=10, export_bin=FALSE, export_track=FALSE) {
	
	
	#Read sam allignment file
	catTime("Reading alignment file [", bam.file, "]",  e={
				what <- c("rname", "strand", "pos",  "qwidth", "mapq")
				flag <- scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = F, 
						hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA, 
						isFirstMateRead = NA, isSecondMateRead = NA, isPrimaryRead = NA, 
						isValidVendorRead = NA, isDuplicate = NA)
				param <- ScanBamParam(flag = flag, simpleCigar = F, what = what)
				
				aln2 <- scanBam(bam.file, param=param)
			})
	
	catTime("Processing alignment file:", e={	
				#Filer ouut reads of quality score lower than "quality_cutoff" [10]
				#aln2 <- (data.frame(aln2)[aln2[[1]]$mapq >= quality_cutoff,]);
				lg <- aln2[[1]]$mapq >= quality_cutoff
				#Construct GRanges object
				#ranges.raw <- GRanges(seqnames = aln2$rname, ranges = IRanges(aln2$pos, width=aln2$qwidth), strand = aln2$strand);
				ranges.raw <- GRanges(seqnames = aln2[[1]]$rname[lg], ranges = IRanges(aln2[[1]]$pos[lg], width=aln2[[1]]$qwidth[lg]), strand = aln2[[1]]$strand[lg]);
				
				#Sort out chromosome naming and sequence lengths (ce6)
				if (genome@provider_version == 'ce6') {	
					seqlevels (ranges.raw) <- seqlevels (genome)[c(1,2,3,4,7,5,6)]
					seqlengths(ranges.raw) <- seqlengths(genome)[c(1,2,3,4,7,5,6)]
				} else if (genome@provider_version == 'hg19') {	
					seqlengths(ranges.raw) <- seqlengths(genome)[seqlevels(ranges.raw)]
				} else if (genome@provider_version == 'mm9') {
					seqlengths(ranges.raw) <- seqlengths(genome)[seqlevels(ranges.raw)]
				} else {
					warning('Unnkonown genome!')
				}
				
				#Resize sequences to 200bp //This resize method can be better with smooth end resizeing (as in peak calling)
				if( !is.null(resize_length) ) { suppressWarnings( ranges.raw <- resize(ranges.raw, resize_length) ) }
				
				#Calculate oryginal BAM file statistics
				num <- countBam(bam.file)
			})
	
	cat(sprintf("\tINFO: %i out of %i [%0.2f%%] reads mapped with %i quality cutoff (%i out of %i [%0.2f%%] nucleotides). \n", length(ranges.raw), num$records, 100*length(ranges.raw)/num$records, quality_cutoff, mean(aln2[[1]]$qwidth[lg])*length(ranges.raw), num$nucleotides, 100*(mean(aln2[[1]]$qwidth[lg])*length(ranges.raw))/num$nucleotides))
	
	rm(aln2)
	
	return(ranges.raw)
}


catTime <- function(..., e=NULL, file="", gc=FALSE) {
	cat(..., "...", sep="", file=file, append=TRUE)
	cat("\t<", system.time(e, gcFirst=gc)[3], "s>\n", sep="", file=file, append=TRUE)	
}