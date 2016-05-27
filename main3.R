library(Rsamtools)
library(GenomicAlignments)

gr <- GRanges("chr14",IRanges(104029449,104029449),"*")
pathToFiles <- "data"
files <- list.files(pathToFiles, full.names=TRUE)

###########
#first version
###########
#part 1 - import to a galignments object
which <- gr  
what <- scanBamWhat()  
flag <- scanBamFlag(isUnmappedQuery = FALSE)
param <- ScanBamParam(flag = flag, which = which, what = what)

bf <- BamFile(files[1])
open(bf)
ga <- readGAlignments(bf, param = param)
close(bf)

#part 2 - summarize allele counts
my_IGPOI <- gr
seqlevels(ga) <- seqlevels(my_IGPOI) 
qseq <- mcols(ga)$seq
			
nuclpiles <- pileLettersAt(qseq, seqnames(ga), start(ga), cigar(ga), my_IGPOI)
nstr <- table(strsplit(as.character(nuclpiles), ""))

nstr


###########
#second version
###########
countF <-
	function(x){
	x[["seq"]][-5,,1]
}                     


which <- gr
flag <- scanBamFlag(isUnmappedQuery = FALSE)

p1 <- ApplyPileupsParam(flag=flag,
						which=which, 
						minBaseQuality = 0L,
						what="seq",
						yieldBy = "position",
						yieldAll=TRUE,
						maxDepth=.Machine$integer.max,
						)

pf <- PileupFiles(bf)
res <- applyPileups(pf, countF, param=p1)



###########
# investigate
###########

#why do these indicate different read depth?
length(ga)
nstr
res


