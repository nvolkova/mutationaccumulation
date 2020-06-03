#####################################################################################
########## script for visualizing complex rearrangements, 14.02.2017 ################
#####################################################################################
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
# ymax - the height of the plot will be max(ymax, max(coverage / coverage in ref))
# unit - bin length for coverage
chr_lens <- c(15072434,15279421,13783801,17493829,20924180,17718942)
names(chr_lens) <- c('I','II','III','IV','V','X')

plot_sv_new <- function(SV,ref.reads=44773320,ymax=5,unit =100, 
                        REFPATH="/PATH/TO/REF/BW", 
                        FILEPATH="/PATH/TO/SAMPLE/BW",
                        ending = '.bw', 
                        main = NULL) # assumes that bw files have the names like "sample_name.bw"
{
  COLOR=c("#E7298A","#66A61E","#E6AB02","cyan")
  names(COLOR) = c("DUP","DEL","INV","BND")
  
  file.ref <- paste(REFPATH) # path to reference BIGWIG
  #file <- paste(FILEPATH,as.character(SV$Sample[1]),"/",as.character(SV$Sample[1]),ending,sep="") # path to BIGWIG of the sample of interest
  file <- paste(FILEPATH, as.character(SV$Sample[1]),ending,sep="") 
  file_bw <- paste("~/bamstats/",as.character(SV$Sample[1]),".stat.dat",sep="") # path to the bamstat file of interest
  
  alt.reads <- sum(read.table(file_bw)$V3[1:6])
  r.ratio <- ref.reads / alt.reads
  
  chr <- as.character(SV$CHR1[1])
  coords <- c(min(as.numeric(SV$POS1)),max(as.numeric(SV$POS2)))
  big.ranges <- c(max(1,coords[1] - (coords[2]-coords[1])),min(coords[2] + (coords[2]-coords[1]),chr_lens[chr]))
  x <- seq(from=big.ranges[1], to=big.ranges[2], by = unit) #???
  big.ranges <- GRanges(seqnames = chr, ranges = IRanges(start=big.ranges[1], end=big.ranges[2]))
  region <- import(BigWigFile(file),which=big.ranges)
  region.ref <- rtracklayer::import(file.ref,which=big.ranges)
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
  
  ymax = 5
  pdf('~/CD0009f.pdf',7,5)
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
  dev.off()
}


# 2 chroms

# SV - table with breakpoints (from svmat)
# ymax - the height of the plot will be max(ymax, max(coverage / coverage in ref))
# unit1 - bin length for coverage in first chromosome
# unit2 - bin length for coverage in second chromosome

plot_sv_int_new <- function(SV,ref.reads=44773320,ymax=5, unit1 = 100, unit2 = 100, 
                            REFPATH="PATH/TO/REFERENCE/BW", 
                            FILEPATH="PATH/TO/SAMPLE/BW",
                            ending = '.bw') # assumes that bw files have the names like "sample_name.bw"
{
  COLOR=c("#E7298A","#66A61E","#E6AB02","cyan")
  names(COLOR) = c("TD","DEL",'INV','BND')
  
  file.ref <- REFPATH # path to reference BW
  file <- paste(FILEPATH,as.character(SV$Sample[1]),ending,sep="") # path to BW of interest
  file_bw <- paste("~/bamstats/",as.character(SV$Sample[1]),".stat.dat",sep="") # path to bamstats
  
  alt.reads <- sum(read.table(file_bw)$V3[1:6])
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
  region1 <- import(file,which=big.ranges1)
  region2 <- import(file,which=big.ranges2)
  region.ref1 <- import(file.ref,which=big.ranges1)
  region.ref2 <- import(file.ref,which=big.ranges2)
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
  
  pdf('CD0751d_variant_4_upd_real.pdf',7,5)
  plot(NA,NA,xlim = c(0, 100), ylim = c(0,ymax), bty = 'n', xaxt = 'n', xlab = 'Position', ylab= 'Coverage',
       main = paste0('Complex variant for ', as.character(unique(SV$Sample)), ', chromosomes ',paste0(chrs, collapse = ',')))    
  axis(side = 1, lwd = 0.3, at = c(0,45,55,100),labels = paste0(round(c(x1[1],x1[length(x1)],x2[1],x2[length(x2)])/1000),' kb'),
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
plot_sv_int(fullsvmat[fullsvmat$NAME == 'CD0009f' & fullsvmat$CLUSTER == 4,], unit1 = 10000, unit2 = 50)
dev.off()

pdf('Complex variant in CD0373i.pdf',7,4)
plot_sv(fullsvmat[fullsvmat$NAME == 'CD0373i' & fullsvmat$CLUSTER == 1,])
dev.off()

pdf('Complex variant in CD0374i(dog-1).pdf',7,4)
plot_sv_int(fullsvmat[fullsvmat$NAME == 'CD0374i' & fullsvmat$CLUSTER == 3,], unit1=100, unit2=10)
dev.off()

pdf('Complex variant in CD0375d(exo-3).pdf',7,4)
plot_sv(fullsvmat[fullsvmat$NAME == 'CD0375d' & fullsvmat$CLUSTER == 2,])
dev.off()
