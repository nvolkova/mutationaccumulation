#####################################################################################
########## script for visualizing complex rearrangements             ################
########## N.Volkova, EMBL-EBI, 2019                                 ################
#####################################################################################

# libraries
library(reshape2)
library(ggplot2)
library(VariantAnnotation)
library(rtracklayer)

#####################################################################################

# load rearrangements and correct worm names

data <- read.csv("Sample_annotation_table.csv") # table with worm genotypes
data$Sample <- as.character(data$Sample)
data$Genotype.new <- as.character(data$Genotype.new)
data$Code <- as.character(data$Code)
CD2Mutant <- sapply(data$Code, function(x) {
  t <- unlist(strsplit(x,split="[:]"))
  t[t=="NA"] <- ""
  if (t[4]!="") # mutagen exposure
    return(paste(t[3],substr(t[4],1,3),t[5],t[7],sep=":")) # genotype, mutagen, dose, experiment type, generation
  if (t[4]=="") # mutation accumulation
    return(paste(t[3],t[7],sep=":")) # genotype, experiment type, generation
})
names(CD2Mutant) <- data$Sample

#####################################################################################

# Single chromosome complex events
# SV - table with breakpoints (from svmat)
# ref.reads - total number of reads in the reference file (default - number of reads in CD0850b)
# ymax - the height of the plot will be max(ymax, max(coverage / coverage in ref))
# unit - bin length for coverage
# REFPATH - path to the bigwig track of the reference sample
# FILEPATH - path to the bigwig track of the sample of interest
# STATPATH - path to the bamstat table for the sample of interest
# main - custom plot title, default is NULL

chr_lens <- c(15072434,15279421,13783801,17493829,20924180,17718942) # fixed the info from C. elegans WB235 (Ensembl)
names(chr_lens) <- c('I','II','III','IV','V','X')

plot_sv <- function(SV,ref.reads=44773320,ymax=5,unit =100, 
                        REFPATH="/PATH/TO/REF/BW", 
                        FILEPATH="/PATH/TO/SAMPLE/BW",
                        STATPATH="/PATH/TO/SAMPLE/BAMSTAT",
                        main = NULL)
{
  COLOR=c("#E7298A","#66A61E","#E6AB02","cyan")
  names(COLOR) = c("DUP","DEL","INV","BND")
  
  alt.reads <- sum(read.table(STATPATH)$V3[1:6])
  r.ratio <- ref.reads / alt.reads
  
  chr <- as.character(SV$CHR1[1])
  coords <- c(min(as.numeric(SV$POS1)),max(as.numeric(SV$POS2)))
  big.ranges <- c(max(1,coords[1] - (coords[2]-coords[1])),min(coords[2] + (coords[2]-coords[1]),chr_lens[chr]))
  x <- seq(from=big.ranges[1], to=big.ranges[2], by = unit) #???
  big.ranges <- GRanges(seqnames = chr, ranges = IRanges(start=big.ranges[1], end=big.ranges[2]))
  region <- import(BigWigFile(FILEPATH),which=big.ranges)
  region.ref <- rtracklayer::import(REFPATH,which=big.ranges)
  bins <- GRanges(seqnames=chr,
                  ranges=IRanges(start=x[-length(x)]+1,end=x[-1],names = rep(chr,length(x)-1)),
                  seqinfo=seqinfo(region))
  numvar <- mcolAsRleList(x=region,varname="score")
  numvar.ref <- mcolAsRleList(x=region.ref,varname="score")
  points <- as.numeric(binnedAverage(bins, numvar, varname="score", na.rm = T)$score) + 0.1
  points.ref <- as.numeric(binnedAverage(bins, numvar.ref, varname="score", na.rm = T)$score) + 0.1
  
  if (is.null(main))
    main = paste("Complex rearrangement in ",as.character(unique(SV$Sample))," (",
                 CD2Mutant[as.character(unique(SV$Sample))],") on chromosome ",
                 as.character(SV$CHR1[1]),sep="")
  
  plot(x[-length(x)], points/points.ref * r.ratio, 
       cex=0.4,lwd=2, 
       ylim = c(0,max(ymax,max(points/points.ref * r.ratio,na.rm = T))),# * r.ratio,na.rm=T))),
       xlab="Genomic region", ylab = "Coverage fold-change relative to normal", 
       main = main, bty = 'n')#, xlim = c(1080000,1108000))
  ymax = max(points/points.ref, na.rm = T)
  lowest <- ymax / nrow(SV)
  for (j in 1:nrow(SV)) {
    height <- lowest + (j-1)*lowest
    color <- COLOR[as.character(SV$TYPE)[j]]
    lines(x = c(SV$POS1[j],SV$POS1[j]),
          y = c(0,height),
          col = color, lwd = 5)
    lines(x = c(SV$POS2[j],SV$POS2[j]),
          y = c(0,height),
          col = color, lwd = 5)
    lines(x=c(SV$POS1[j],SV$POS2[j]),
          y=rep(height,2), col=color, lwd=3)
  }
  legend('topleft',cex = 0.5, legend = names(COLOR), fill = COLOR)

}

# 2-chromosome complex events
# SV - table with breakpoints (from svmat)
# ref.reads - total number of reads in the reference file (default - number of reads in CD0850b)
# ymax - the height of the plot will be max(ymax, max(coverage / coverage in ref))
# unit1 - bin length for coverage in first chromosome
# unit2 - bin length for coverage in second chromosome
# REFPATH - path to the bigwig track of the reference sample
# FILEPATH - path to the bigwig track of the sample of interest
# STATPATH - path to the bamstat table for the sample of interest
# main - custom plot title, default is NULL

plot_sv_int <- function(SV,ref.reads=44773320,ymax=5, unit1 = 100, unit2 = 100, 
                            REFPATH="PATH/TO/REFERENCE/BW", 
                            FILEPATH="PATH/TO/SAMPLE/BW",
                            STATPATH="/PATH/TO/SAMPLE/BAMSTAT",
                            ending = '.bw') 
{
  COLOR=c("#E7298A","#66A61E","#E6AB02","cyan")
  names(COLOR) = c("TD","DEL",'INV','BND')
  
  alt.reads <- sum(read.table(STATPATH)$V3[1:6])
  r.ratio <- ref.reads / alt.reads
  
  chrs <- unique(c(as.character(SV$CHR1),as.character(SV$CHR2)))
  coords1 <- c(min(c(as.numeric(SV$POS1[SV$CHR1 == chrs[1]]),as.numeric(SV$POS2[SV$CHR2 == chrs[1]]))),
               max(c(as.numeric(SV$POS1[SV$CHR1 == chrs[1]]),as.numeric(SV$POS2[SV$CHR2 == chrs[1]]))))
  coords2 <- c(min(c(as.numeric(SV$POS1[SV$CHR1 == chrs[2]]),as.numeric(SV$POS2[SV$CHR2 == chrs[2]]))),
               max(c(as.numeric(SV$POS1[SV$CHR1 == chrs[2]]),as.numeric(SV$POS2[SV$CHR2 == chrs[2]]))))
  big.ranges1 <- c(max(1,coords1[1] - 500),min(coords1[2] + 500,chr_lens[chrs[1]]))
  if (coords1[2] - coords1[1] == 0) 
    big.ranges1 <- c(max(1,coords1[1] - 2000),min(coords1[2] + 2000,chr_lens[chrs[1]]))
  big.ranges2 <- c(max(1,coords2[1] - (coords2[2]-coords2[1])/4),min(coords2[2] + (coords2[2]-coords2[1])/4,chr_lens[chrs[2]]))
  if (coords2[2] - coords2[1] == 0) 
    big.ranges2 <- c(max(1,coords2[1] - 2000),min(coords2[2] + 2000,chr_lens[chrs[2]]))
  x1 <- seq(from=big.ranges1[1], to=big.ranges1[2], by = unit1)
  x2 <- seq(from=big.ranges2[1], to=big.ranges2[2], by = unit2) 
  big.ranges1 <- GRanges(seqnames = chrs[1], ranges = IRanges(start=big.ranges1[1], end=big.ranges1[2]))
  big.ranges2 <- GRanges(seqnames = chrs[2], ranges = IRanges(start=big.ranges2[1], end=big.ranges2[2]))
  region1 <- import(FILEPATH,which=big.ranges1)
  region2 <- import(FILEPATH,which=big.ranges2)
  region.ref1 <- import(REFPATH,which=big.ranges1)
  region.ref2 <- import(REFPATH,which=big.ranges2)
  bins1 <- GRanges(seqnames=chrs[1],
                   ranges=IRanges(start=x1[-length(x1)]+1,end=x1[-1],names = rep(chrs[1],length(x1)-1)),
                   seqinfo=seqinfo(region1))
  bins2 <- GRanges(seqnames=chrs[2],
                   ranges=IRanges(start=x2[-length(x2)]+1,end=x2[-1],names = rep(chrs[2],length(x2)-1)),
                   seqinfo=seqinfo(region2))
  numvar1 <- mcolAsRleList(x=region1,varname="score")
  numvar.ref1 <- mcolAsRleList(x=region.ref1,varname="score")
  points1 <- as.numeric(binnedAverage(bins1,numvar1,varname="score",na.rm = T)$score) + 0.1
  points.ref1 <- as.numeric(binnedAverage(bins1,numvar.ref1,varname="score",na.rm = T)$score) + 0.1
  numvar2 <- mcolAsRleList(x=region2,varname="score")
  numvar.ref2 <- mcolAsRleList(x=region.ref2,varname="score")
  points2 <- as.numeric(binnedAverage(bins2,numvar2,varname="score",na.rm = T)$score) + 0.1
  points.ref2 <- as.numeric(binnedAverage(bins2,numvar.ref2,varname="score",na.rm = T)$score) + 0.1
  
  if (is.null(main))
    main = paste0('Complex variant for ', as.character(unique(SV$Sample)), 
                  ', chromosomes ',paste0(chrs, collapse = ','))
  
  plot(NA,NA,xlim = c(0, 100), ylim = c(0,ymax), bty = 'n', 
       xaxt = 'n', xlab = 'Position', ylab= 'Coverage', main = main)   
  axis(side = 1, lwd = 0.3, at = c(0,45,55,100),
       labels = paste0(round(c(x1[1],x1[length(x1)],x2[1],x2[length(x2)])/1000),' kb'),
       las = 2, cex.axis = 0.4)
  abline(v = 50, lty = 2)
  coords1 <- seq(0,45,length.out = length(x1)-1)
  coords2 <- seq(55,100,length.out = length(x2)-1)
  points(coords1, points1/points.ref1 * r.ratio, cex=0.4,lwd=2)
  points(coords2, points2/points.ref2 * r.ratio, cex=0.4,lwd=2)
  text(x = c(25,75), y = rep(ymax,2), labels = chrs)
  for (j in 1:nrow(SV)) {
    color <- COLOR[as.character(SV$TYPE)[j]]
    if (SV$CHR1[j] == chrs[1]) {
      height <- (ymax/2) + (j-1) * (ymax/4)#mean((points1/points.ref1 * r.ratio)[which(abs(x1-SV$POS1[j])<unit1)[1]:which(abs(x1-SV$POS1[j])<unit1)[2]]) + j
      lines(x = rep(coords1[which(abs(x1-SV$POS1[j])<unit1)[1]],2), y = c(0,height), col = color, lwd = 5)
      start.p = coords1[which(abs(x1-SV$POS1[j])<unit1)[1]]
    } else {
      height <- mean((points2/points.ref2 * r.ratio)[which(abs(x2-SV$POS1[j])<unit2)[1]:which(abs(x2-SV$POS1[j])<unit2)[length(which(abs(x2-SV$POS1[j])<unit2))]]) + j
      lines(x = rep(coords2[which(abs(x2-SV$POS1[j])<unit1)[1]],2), y = c(0,height), col = color, lwd = 5)
      start.p = coords2[which(abs(x2-SV$POS1[j])<unit2)[1]]
    }
    if (SV$CHR2[j] == chrs[1]) {
      #height <- mean((points1/points.ref1 * r.ratio)[which(abs(x1-SV$POS2.1[j])<unit1)[1]:which(abs(x1-SV$POS2.2[j])<unit1)[1]])
      lines(x = rep(coords1[which(abs(x1-SV$POS2[j])<unit1)[1]],2), y = c(0,height), col = color, lwd = 5)
      end.p = coords1[which(abs(x1-SV$POS2[j])<unit1)[1]]
    } else {
      #height <- mean((points2/points.ref2 * r.ratio)[which(abs(x2-SV$POS2.1[j])<unit2)[1]:which(abs(x2-SV$POS2.2[j])<unit2)[1]])
      lines(x = rep(coords2[which(abs(x2-SV$POS2[j])<unit1)[1]],2), y = c(0,height), col = color, lwd = 5)
      end.p = coords2[which(abs(x2-SV$POS2[j])<unit2)]#[2]]
    }
    lines(x = c(start.p, end.p), y=rep(height,2),col=color, lwd=3)
  }
  legend('topleft', legend = names(COLOR), fill = COLOR)
  dev.off()
}


#########
# Usage examples - sapply a table of breakpoints for a single COMPLEX event
# svmat should look e.g. like this:
#             CHR1     POS1 CHR2     POS2 READS TYPE  Sample clust.type
# INV00000953    V  6332277    V  6332849    36  INV CD0394d    COMPLEX
# INV00000954    V  6334532    V  6335418    11  INV CD0394d    COMPLEX
# BND00000955    V  6336793    I 10124012    21  BND CD0394d    COMPLEX
# BND00000956    V  6337422    I 10030202    12  BND CD0394d    COMPLEX
# BND00000957    V  6338031    I 10032498    15  BND CD0394d    COMPLEX

pdf('Complex variant in CD0009f.pdf',7,4)
plot_sv_int(fullsvmat[fullsvmat$NAME == 'CD0009f' & fullsvmat$CLUSTER == 4,], 
            unit1 = 10000, unit2 = 50,
            REFPATH='/path/to/CD0001b.bw',
            FILEPATH='/path/to/CD0009f.bw',
            STATPATH='/path/to/CD0009f.bamstat.txt')
dev.off()


pdf('Complex variant in CD0375d(exo-3).pdf',7,4)
plot_sv(fullsvmat[fullsvmat$NAME == 'CD0375d' & fullsvmat$CLUSTER == 2,],
        REFPATH='/path/to/CD0001b.bw',
        FILEPATH='/path/to/CD0375d.bw',
        STATPATH='/path/to/CD0375d.bamstat.txt')
dev.off()

# sessionInfo()
# R version 3.6.3 (2020-02-29)
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] rtracklayer_1.44.4          VariantAnnotation_1.30.1    Rsamtools_2.0.3             Biostrings_2.52.0          
# [5] XVector_0.24.0              SummarizedExperiment_1.14.1 DelayedArray_0.10.0         BiocParallel_1.18.1        
# [9] matrixStats_0.55.0          Biobase_2.44.0              GenomicRanges_1.36.1        GenomeInfoDb_1.20.0        
# [13] IRanges_2.18.3              S4Vectors_0.22.1            BiocGenerics_0.30.0         ggplot2_3.2.1              
# [17] reshape2_1.4.3             
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.2               lattice_0.20-38          prettyunits_1.0.2        assertthat_0.2.1         digest_0.6.22           
# [6] R6_2.4.0                 plyr_1.8.4               RSQLite_2.1.2            evaluate_0.14            httr_1.4.1              
# [11] pillar_1.4.2             zlibbioc_1.30.0          rlang_0.4.5              GenomicFeatures_1.36.4   progress_1.2.2          
# [16] lazyeval_0.2.2           rstudioapi_0.10          blob_1.2.0               Matrix_1.2-18            rmarkdown_1.16          
# [21] stringr_1.4.0            RCurl_1.95-4.12          bit_1.1-14               biomaRt_2.40.5           munsell_0.5.0           
# [26] compiler_3.6.3           xfun_0.10                pkgconfig_2.0.3          htmltools_0.4.0          tidyselect_0.2.5        
# [31] tibble_2.1.3             GenomeInfoDbData_1.2.1   XML_3.98-1.20            crayon_1.3.4             dplyr_0.8.3             
# [36] withr_2.1.2              GenomicAlignments_1.20.1 bitops_1.0-6             grid_3.6.3               gtable_0.3.0            
# [41] DBI_1.0.0                magrittr_1.5             scales_1.0.0             stringi_1.4.3            vctrs_0.2.4             
# [46] tools_3.6.3              bit64_0.9-7              BSgenome_1.52.0          glue_1.3.1               purrr_0.3.3             
# [51] hms_0.5.3                yaml_2.2.0               AnnotationDbi_1.46.1     colorspace_1.4-1         memoise_1.1.0           
# [56] knitr_1.25   
#
