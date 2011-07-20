# TODO: Add comment
# 
# Author: przemol
###############################################################################
#setwd('/Volumes/raid0/rBeadsNew/Tests')

#require(compiler)
#enableJIT(level=0)

#require(digest)
#reName <- function(file_name, proccesing="norm", scale='linear', resolution='25bp', ext='.wig', prefix='C', uid=1) {
#	if ( grepl("^([A-Za-z0-9]+)\\^([A-Za-z0-9]+)_([A-Za-z0-9\\-]+)\\^([A-Za-z0-9]+)\\^([A-Za-z0-9]+)\\^([A-Za-z0-9]+)_([A-Za-z0-9]+)\\^([A-Za-z0-9]+)\\^([A-Za-z0-9]+)", file_name, perl=T) ) {
#		
#		ff <- strsplit( unlist(strsplit(file_name, '\\.', perl=T)) , '(\\^|_)', perl=T)[[1]]
#		out <- sprintf('%s^%s_%s^%s^%s^%s_%s^%s^%s_%s', ff[1], ff[2], ff[3], ff[4], ff[5], ff[6], proccesing, scale, resolution, ff[10])
#		md5 <- substr(digest(out, algo='md5'), 1, 2)
#		num = sprintf("%05d", uid)
#		return(sprintf('%s^%s%s%s%s', out, prefix, md5, num, ext))
#	} else {
#		warning('Wrong name format')
#		return(sprintf("%s_%s^%s^%s%s", unlist(strsplit(file_name, "\\."))[1], proccesing, scale, resolution, ext))
#	}
#}
#ParseName <- function(file_name, proccesing="norm", file_type='wiggle', ext='.wig', prefix='C', uid=1) {
#	if ( grepl("^([A-Za-z0-9]+)\\^([A-Za-z0-9]+)_([A-Za-z0-9\\-]+)\\^([A-Za-z0-9]+)\\^([A-Za-z0-9]+)\\^([A-Za-z0-9]+)_([A-Za-z0-9]+)\\^([A-Za-z0-9]+)\\^([A-Za-z0-9]+)", file_name, perl=T) ) {
#		
#		ff <- strsplit( unlist(strsplit(file_name, '\\.', perl=T)) , '(\\^|_)', perl=T)[[1]]
#		ext <- strsplit( unlist(strsplit(file_name, '\\.', perl=T)) , '(\\^|_)', perl=T)[[2]]
#		df <- data.frame(Factor=I(ff[1]), Antibody=I(ff[2]), ExtractID=I(ff[3]), Crosslinker=I(ff[4]), Strain=I(ff[5]), Stage=I(ff[6]), Processing=I(ff[7]), Resolution=I(ff[8]), Scale=I(ff[9]), ContactExpID=I(ff[10]), UID=I(ff[11]), Extension=I(ext))
#		return(df)
#	} else {
#		warning('Wrong name format')
#	}
#}

sourceDir <- function(path, trace = TRUE, ...) {
	for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
		if(trace) cat(nm,":")           
		source(file.path(path, nm), ...)
		if(trace) cat("\n")
	}
}
source.dir <- '/Users/przemol/git/rBEADS_v2.0/rBEADS_v2.0/source'
sourceDir(source.dir)
#load(file.path(source.dir, 'precalculated/MappabilityCe6_v2.Rdata'))
#load(file.path(source.dir, 'precalculated/nonMappableFilterCe6_v1.Rdata'))
#load(file.path(source.dir, 'precalculated/SummedFormaldehydeInput_v2.Rdata'))
#load(file.path(source.dir, 'precalculated/SummedEGSInput_v2.Rdata'))
rm(source.dir)

require(BSgenome.Celegans.UCSC.ce6)
require(BSgenome.Hsapiens.UCSC.hg19)
require(BSgenome.Mmusculus.UCSC.mm9)
require(BSgenome.Dmelanogaster.UCSC.dm3)


beads <- function(config='BeadsConfig.csv') {
	
	data.dir <- getwd()
	setwd(data.dir)
	files <- suppressWarnings(read.csv(config))
	
	genome <- get(as.character(files$genome))
	MAP <- load(as.character(files$maptrack))
	MAPF <- load(as.character(files$mapfilter))
	ext <- as.integer(as.character(files$extlength))
	GCcutoff <- c(as.integer(files$gclow), as.integer(files$gchigh))
	
	#MAP <- data('MappabilityCe6', package='rBEADS')
	#MAPF <- data('nonMappableFilterCe6', package='rBEADS')
	#FRMI <- load(file.path(source.dir, 'precalculated/SummedFormaldehydeInput_v2.Rdata'))

	for(i in 1:nrow(files))  {
		
		setwd(data.dir)
		
		##Set up naming
		sample.d <- reName(as.character(files$Sample[i]), 'dir', 'na', 'na', '', 'D', 0)
		
		
		##Log output
		dir.create(sample.d); setwd(sample.d)
		sink(file = 'LOG.txt', type = "output", split=TRUE)
		##zz <- file("LOG.ERROR.txt", open="wt"); sink(zz, split=TRUE); sink(zz, type="message")
		
		cat("/* File created on", date(), "*/\n")
		cat(i, ") Processing: ", as.character(files$Sample[i]), ' & ', as.character(files$Control[i]), "\n", sep='')
		
		## if(files$Control[i] == 'FRM') {
		##     cat('INFO: Using summed formaldehyde input!\n')
		##     if(exists('summed.frm')) { control.map <- summed.frm } else { control.map <- get(data('summed.frm', package='rBEADS')) }
		## } else if(files$Control[i] == 'EGS') {
		##     cat('INFO: Using summed EGS input!\n')
		##     if(exists('summed.egs')) { control.map <- summed.egs } else { control.map <- get(data('summed.egs', package='rBEADS')) }
		## }  else if(grepl('Rdata$', as.character(files$Control[1]))) {
		##     cat('INFO: Using summed', as.character(files$Control[1]), 'input!\n')
		##     control.map <- get(load( file.path(data.dir, as.character(files$Control[1])) ))
		## } else {
			##Control [INPUT]
			control.d <- reName(as.character(files$Control[i]), 'na', 'na', 'na', '', 'D', 0)
			if (grepl('bam$', as.character(files$Control[1]))) {
				#Valid BAM control for all genomes!
				control.re <- ImportBAM(bam.file=file.path(data.dir, as.character(files$Control[i])), genome=genome, desc=control.d, resize_length=ext, quality_cutoff=10)
			} else if (grepl('export.fq', as.character(files$Control[1]))) {
				control.re <- ImportSolexa(file=as.character(files$Control[i]), data.dir=data.dir, desc=control.d, resize_length=200, quality_cutoff=10, export_bin=FALSE, export_track=FALSE)
			} else if (grepl('Rdata$', as.character(files$Control[1]))) {
				sample.re <- get(load( file.path(data.dir, as.character(files$Control[i])) ))
			} else {
				stop('Unknown input file type for sample!')
			}
			control.gc <- GCCorrection(ranges.raw=control.re, enriched_regions=NULL, nonMappableFilter=get(MAPF), genome=genome, desc=control.d, smoothing_spline=FALSE, cutoff=GCcutoff)
			control.map <- MappabilityCorrection(GCnormTrack=control.gc, mappabilityTrack=get(MAP), genome=genome)
		## } 
		
		##Sample [ChIP]
		if (grepl('bam$', as.character(files$Sample[1]))) {
			#Valid BAM sample for all genomes!
			sample.re <- ImportBAM(bam.file=file.path(data.dir, as.character(files$Sample[i])), genome=genome, desc=sample.d, resize_length=ext, quality_cutoff=10)
		} else if (grepl('export.fq', as.character(files$Sample[1]))) {
			sample.re <- ImportSolexa(file=as.character(files$Sample[i]), data.dir=data.dir, desc=sample.d, resize_length=200, quality_cutoff=10, export_bin=FALSE, export_track=FALSE)
		} else if (grepl('Rdata$', as.character(files$Sample[1]))) {
			sample.re <- get(load( file.path(data.dir, as.character(files$Sample[i])) ))
		} else {
			stop('Unknown input file type for sample!')
		}
		if(is.null(files$ER[i]) |  files$ER[i]=='auto' ){ 
			sample.er <- EnrichedRegions(sample.re, desc=sample.d, genome=genome) 
		} else if (files$ER[i]=='no') { 
			cat('INFO: Skipping ER auto finder!\n')
			sample.er <- NULL
		} else {
			sample.er <- import.bed(files$ER[i], genome = "hg18")
		}
		sample.gc <- GCCorrection(ranges.raw=sample.re, enriched_regions=sample.er, nonMappableFilter=get(MAPF), genome=genome, desc=sample.d, smoothing_spline=FALSE, cutoff=GCcutoff)
		sample.map <- MappabilityCorrection(GCnormTrack=sample.gc, mappabilityTrack=get(MAP), genome=genome)
		
		##Normalization
		sample.norm <- DivStep(sample.map, sample.map, sample.er, genome=genmome)
		
		
		
		#EXPORT	
		dir.create('Rbinaries'); setwd('Rbinaries')
			sample.re.cov <- coverage(sample.re)
			save(sample.re.cov, file=reName(sample.d, 'Raw', 'linar', '01bp', ext='.Rdata', 'B', 0))
			save(sample.gc, 	file=reName(sample.d, 'GCCorrected', 'linar', '01bp', ext='.Rdata', 'B', 0))
			save(sample.map, 	file=reName(sample.d, 'MapCorrected', 'linar', '01bp', ext='.Rdata', 'B', 0))
			save(sample.norm, 	file=reName(sample.d, 'NORM', 'linar', '01bp', ext='.Rdata', 'B', 0))
		setwd('..')
		
		if (exists('sample.re')) {
			dir.create('ChIP_QC'); setwd('ChIP_QC')
				BinTrack(coverage(sample.re), n=25, smooth=FALSE, out=reName(sample.d, 'Raw', 'linear', '25bp', ext='.wig'), type="WIG", name='RAW alignment', col='darkred')
				BinTrack(sample.gc, n=25, smooth=FALSE, out=reName(sample.d, 'GCCorrected', 'linear', '25bp', ext='.wig'), type="WIG", name='GC corrected', col='darkblue')
				BinTrack(sample.map, n=25, smooth=FALSE, out=reName(sample.d, 'MapCorrected', 'linear', '25bp', ext= '.wig'), type="WIG", name='Mappability corrected', col='blue')
			setwd('..')
		}
		if (exists('control.re')) {
			dir.create('Input_QC'); setwd('Input_QC')
				BinTrack(coverage(control.re), n=25, smooth=FALSE, out=reName(control.d, 'Raw', 'linear', '25bp', ext='.wig'), type="WIG", name='RAW alignment', col='black')
				BinTrack(control.gc, n=25, smooth=FALSE, out=reName(control.d, 'GCCorrected', 'linear', '25bp', ext='.wig'), type="WIG", name='GC corrected', col='darkblue')
				BinTrack(control.map, n=25, smooth=FALSE, out=reName(control.d, 'MapCorrected', 'linear', '25bp', ext='.wig'), type="WIG", name='Mappability corrected', col='blue')
			setwd('..')
		}
		dir.create('NORMALIZED'); setwd('NORMALIZED')
			BinTrack(sample.norm, n=25, smooth=FALSE, out=reName(sample.d, 'NORM', 'linear', '25bp', ext='.wig', 'N', 0), type="WIG", name=paste(sample.d, 'BEADS Normalised'), col='darkgreen')
			BinTrack(sample.norm, n=25, smooth=FALSE, out=reName(sample.d, 'NORM', 'Log2', '25bp', ext='.wig', 'N', 0), type="WIG", zscore=FALSE, log2=TRUE, name=paste(sample.d, 'BEADS Normalised & log2'), col='darkgreen')
			BinTrack(sample.norm, n=25, smooth=FALSE, out=reName(sample.d, 'NORM', 'Zscore', '25bp', ext='.wig', 'N', 0), type="WIG", zscore=TRUE, log2=FALSE, name=paste(sample.d, 'BEADS Normalised & z-scored'), col='darkgreen')
			BinTrack(sample.norm, n=25, smooth=FALSE, out=reName(sample.d, 'NORM', 'linear', '25bp', ext='.gff', 'N', 0), type="GFF")
		setwd('..')
				
		cat("/* Calculation complited on", date(), "*/\n")
		##sink(type="message"); sink(); close(zz); 
		sink(); setwd('..');
	}
}

if (0) {
	setwd('/Volumes/raid0/rBeadsNew/Rpackage')
	
	system('rm -R rBEADS', intern=TRUE)
	package.skeleton('rBEADS')
	
	file.remove('rBEADS/DESCRIPTION')
	file.copy('DESCRIPTION', 'rBEADS/')
	file.remove(dir('rBEADS/man/', full.name=TRUE))
	file.copy(dir('man/', full.name=TRUE), 'rBEADS/man/')
	system('R CMD build rBEADS', intern=TRUE)
	
	file.remove('/Library/WebServer/Documents/rBEADS/src/contrib/PACKAGES')
	file.copy('DESCRIPTION', '/Library/WebServer/Documents/rBEADS/src/contrib/PACKAGES')
	file.remove('/Library/WebServer/Documents/rBEADS/src/contrib/rBEADS_1.2.2.tar.gz')
	##file.remove('/Library/WebServer/Documents/rBEADS/rBEADS_1.2.1.tgz')
	file.copy('rBEADS_1.2.2.tar.gz', '/Library/WebServer/Documents/rBEADS/src/contrib/')
	##file.copy('rBEADS_1.2.1.tar.gz', '/Library/WebServer/Documents/rBEADS/rBEADS_1.2.1.tgz')
	
	#Direct install
	remove.packages('rBEADS')
	#install.packages('rBEADS', contriburl='http://ja-mac1.gurdon.cam.ac.uk/rBEADS', type='source')
	#R CMD Rd2pdf man/beads.Rd --force
	install.packages('rBEADS', repos=c(rBEADS='http://ja-mac1.gurdon.cam.ac.uk/rBEADS', rforge="http://r-forge.r-project.org", CRAN="http://www.stats.bris.ac.uk/R"), type='source', dependencies=TRUE)
	
	##REPOS install
	source("http://bioconductor.org/biocLite.R")
	biocLite(c('Rsamtools', 'BSgenome.Celegans.UCSC.ce6', 'digest', 'rtracklayer', 'GenomicRanges'))
	install.packages('rBEADS', repos=c(rBEADS="http://ja-mac1.gurdon.cam.ac.uk/rBEADS", rforge="http://r-forge.r-project.org", getOption("repos")), type='source', dependencies=TRUE)
}
#
#
#terz dodam jakies bzdury na koncu pliku

	
