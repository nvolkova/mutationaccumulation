################################################################################
## Collecting and visualising information about SVs in different C. elegans   ##
## mutation accumulation lines                                                ##
## N.Volkova, EMBL-EBI, 2019                                                  ##
################################################################################

# libraries
library(ggplot2)
library(reshape2)
library(rtracklayer)
library(BSgenome)
worm_ref_genome <- "BSgenome.Celegans.UCSC.ce11"
library(worm_ref_genome, character.only = TRUE)
source('useful_functions.R')
library(vioplot)
library(RColorBrewer)

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

# load the table with SVs, available upon request
SVclust.new <- list()
PATH_TO_SV_TABLES='/path/to/per/sample/tables/with/SVs'
for (samplename in names(CD2Mutant)) {
  
  SVclust.new <- read.table(paste0(PATH_TO_SV_TABLES, '/',samplename,'.tsv'),
                            sep = '\t', header = T)
  
}

# Select genotypes of interest (those with a significant amount of SVs)
sv.properties <- data.frame(genotype = c('N2', 'atm-1', 'brc-1', 'brc-1,ced-3', 
                                         'brc-1,cep-1', 'dog-1', 'helq-1', 'him-6', 
                                         'him-6,ced-3', 'him-6,cep-1', 'mus-81', 'parp-2', 
                                         'rip-1', 'rtel-1', 'slx-1', 'smc-5', 'smc-6', 'wrn-1'),
                            total_number = rep(0,18),
                            deletions = rep(0,18),
                            tandem_duplications = rep(0,18),
                            svs_with_Grich = rep(0,18),
                            indels_in_Grich = rep(0,18),
                            total_indels = rep(0,18),
                            svs_rep_left = rep(0,18),
                            svs_rep_right = rep(0,18),
                            variants_in_telomeres = rep(0,18),
                            variants_in_repeats = rep(0,18),
                            no.samples = rep(0,18))
sv.properties <- sv.properties[-c(10,12,18),] # genotypes excluded later
sv.properties$genotype <- as.character(sv.properties$genotype)

########################################################################################

# GC rich regions
# from Marsico et al. 2019

GCPATH='/path/to/GC/tracks/'
gc_plus <- import(paste0(GCPATH,'Celegans_all_w15_th-1_plus.hits.max.PDS.w50.35.bed'))
gc_minus <- import(paste0(GCPATH, 'Celegans_all_w15_th-1_minus.hits.max.PDS.w50.35.bed'))
seqlevels(gc_plus) <- sapply(seqlevels(gc_plus), function(x) unlist(strsplit(x,split = '[_]'))[2])
seqlevels(gc_minus) <- sapply(seqlevels(gc_minus), function(x) unlist(strsplit(x,split = '[_]'))[2])

########################################################################################

# Replication directions
# from Pourkarimi et al. 2016

REPPATH='/path/to/reptime/tracks/'

repTimeN2LPRplus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417785_N2L_Pr_plus.wig.gz"), format="WIG")
repTimeN2LPRminus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417785_N2L_Pr_neg.wig.gz"), format="WIG")
repTimeN2LRRplus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417786_N2L_RR_plus.wig.gz"), format="WIG")
repTimeN2LRRminus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417786_N2L_RR_neg.wig.gz"), format="WIG")

repTimeN2PRplus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417781_N2_Pr_plus.wig.gz"), format="WIG")
repTimeN2PRminus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417781_N2_Pr_neg.wig.gz"), format="WIG")
repTimeN2RRplus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417782_N2_RR_plus.wig.gz"), format="WIG")
repTimeN2RRminus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417782_N2_RR_neg.wig.gz"), format="WIG")

repTimeegl30PRplus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417783_egl30_Pr_plus.wig.gz"), format="WIG")
repTimeegl30PRminus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417783_egl30_Pr_neg.wig.gz"), format="WIG")
repTimeegl30RRplus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417784_egl30_RR_plus.wig.gz"), format="WIG")
repTimeegl30RRminus <- import(paste0(REPPATH, '/', "GSE90939_RAW/GSM2417784_egl30_RR_neg.wig.gz"), format="WIG")


# Replication fork directionality
repTime <- granges(repTimeN2LPRplus)
repTime <- resize(repTime, width(repTime) + 99, fix="end")
mcols(repTime) <- DataFrame(N2.1 = (-repTimeN2PRminus$score-repTimeN2PRplus$score)/(-repTimeN2PRminus$score+repTimeN2PRplus$score),
                            N2.2 = (-repTimeN2RRminus$score-repTimeN2RRplus$score)/(-repTimeN2RRminus$score+repTimeN2RRplus$score),
                            elg30.1 =  (-repTimeegl30PRminus$score-repTimeegl30PRplus$score)/(-repTimeegl30PRminus$score+repTimeegl30PRplus$score),
                            elg30.2 =  (-repTimeegl30RRminus$score-repTimeegl30RRplus$score)/(-repTimeegl30RRminus$score+repTimeegl30RRplus$score),
                            N2L.1 =  (-repTimeN2LPRminus$score-repTimeN2LPRplus$score)/(-repTimeN2LPRminus$score+repTimeN2LPRplus$score),
                            N2L.2 =  (-repTimeN2LRRminus$score-repTimeN2LRRplus$score)/(-repTimeN2LRRminus$score+repTimeN2LRRplus$score))
t <- as.matrix(mcols(repTime))
t[is.nan(t)] <- NA
m <- rowMeans(t, na.rm = t)
s <- rowMeans(sqrt((t-m)^2)) # standard errors across different worms
table(abs(m)-2*s > 0)
# FALSE   TRUE 
# 507790 455057 
# 45% of the genome assigned with directionality

repStrand <- granges(repTime)
strand(repStrand)[which(m<0)] <- "-" # leftward moving fork
strand(repStrand)[which(m>0)] <- "+" # rightward moving fork
strand(repStrand)[which(abs(m)-2*s < 0)] <- "*" # equally possible directions
seqlevels(repStrand) <- sub("chr","", seqlevels(repStrand))

########################################################################################

# Repetitive regions
# from www.girinst.org/downloads/repeatmaps/C.Elegans

REPEATS='/path/to/downloaded/repetitive/region/tracks/'
repeat.table <- do.call('rbind', lapply(paste0(REPEATS,'/',
                                               c('CEmapI.gz',
                                                 'CEmapII.gz',
                                                 'CEmapIII.gz',
                                                 'CEmapIV.gz',
                                                 'CEmapV.gz',
                                                 'CEmapX.gz')),read.table,header=F))
head(repeat.table)
colnames(repeat.table) <- c('Chromosome', 'Start_on_chrom', 'End_on_chrom', 
                            'Name_of_similar_in_Repbase', 'Start_RepBase','End_RepBase',
                            'Orientation','Identity')
repeat.table.ranges <- GRanges(seqnames = paste0('chr',substr(as.character(repeat.table$Chromosome),
                                                        4,
                                                        nchar(as.character(repeat.table$Chromosome)))),
                               ranges = IRanges(start = repeat.table$Start_on_chrom, 
                                                end = repeat.table$End_on_chrom, 
                                                names = paste0(repeat.table$Chromosome,':',
                                                               repeat.table$Start_on_chrom,'-',
                                                               repeat.table$End_on_chrom)),
                               seqinfo = seqinfo(get(worm_ref_genome)),
                               strand = c('+','-')[as.numeric(repeat.table$Orientation == 'd') + 1])

mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
seqlevels(repeat.table.ranges) <- substr(as.character(seqlevels(repeat.table.ranges)),4,
                                        nchar(as.character(seqlevels(repeat.table.ranges))))

########################################################################################

# Subtelomeric regions
WBcel235 <- getSeq(BSgenome.Celegans.UCSC.ce11)
chr_sizes <- width(WBcel235)
names(chr_sizes) <- c("I","II","III","IV","V","X","MtDNA")
genome_size = sum(as.numeric(chr_sizes))
# consider 200 bp on each side telomeric, and up to 16kbp from the ends - subtelomeric
subtelomere <- matrix(1,nrow=6,ncol=4,dimnames=list(c("I","II","III","IV","V","X"),c("l.start","l.end","r.start","r.end")))
subtelomere[,1] <- rep(200,1)
subtelomere[,2] <- rep(16000,1)
subtelomere[,3] <- chr_sizes[-7] - 16000
subtelomere[,4] <- chr_sizes[-7] - 200

########################################################################################

# Go through all samples
endings <- c(":5",":6",":10",":20",":15",":40") # different generations

replication_tables_td <- list()
replication_tables_del <- list()

td_size <- list()
del_size <- list()

# Can try but rather do gene-by-gene
for (gene in sv.properties$genotype) {
  
  sv <- do.call('rbind',SVclust.new[names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]])
  no.sv <- do.call('rbind', delly.filtcounts[names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]])
  sv <- sv[sv$clust.type!='some',]
  # unique variants only, duplicated.breakpoints in useful_functions.R
  if (gene == 'N2') {
    sv <- sv[-c(6:7),]
    sv[5,'POS2'] <- 8379643
  }
  if (gene == 'atm-1') {
    sv[4,'POS2'] <- sv[5,'POS2']
    sv[11,'POS2'] <- sv[12,'POS2']
    sv <- sv[-c(5,7,12),]
  }
  if (gene %in% c('brc-1',"brc-1,cep-1", 'rtel-1','helq-1')) {
    sv <- sv[!duplicated.breakpoints(sv),]
  }
  if (gene == "brc-1,ced-3") {
    sv <- sv[-1,]
  }
  if (gene == 'dog-1') {
    sv[8,'POS2'] <- sv[9,'POS2']
    sv[12,'POS2'] <- sv[13,'POS2']
    sv <- sv[-c(9,13),]
    sv <- sv[!duplicated.breakpoints(sv),]
  }
  if (gene == 'him-6') {
    sv[5,'POS2'] <- sv[6,'POS2']
    sv <- sv[-c(6:7,14),]
    sv <- sv[!duplicated.breakpoints(sv),]
  }
  if (gene == 'him-6,ced-3') {
    sv[4,'POS2'] <- sv[5,'POS2']
    sv[6,'POS2'] <- sv[7,'POS2']
    sv[10,'POS2'] <- sv[11,'POS2']
    sv <- sv[-c(5,7,11),]
  }
  if (gene == 'mus-81') {
    sv[16,'POS2'] <- sv[18,'POS2']
    sv[34,'POS2'] <- sv[36,'POS2']
    sv <- sv[-c(3:4,17:18,25:26,35:36,39:40,51:52),]
    sv <- sv[!duplicated.breakpoints(sv),]
  }
  if (gene == 'rip-1') {
    sv[1,'POS2'] <- sv[2,'POS2']
    sv[5,'POS2'] <- sv[6,'POS2']
    sv[7,'POS2'] <- sv[8,'POS2']
    sv[10,'POS2'] <- sv[11,'POS2']
    sv <- sv[-c(2,6,8,11,14,17),]
    sv <- sv[!duplicated.breakpoints(sv),]
  }
  if (gene == 'slx-1') {
    sv[15,'POS2'] <- sv[17,'POS2']
    sv[19,'POS2'] <- sv[20,'POS2']
    sv[37,'POS2'] <- sv[39,'POS2']
    sv[40,'POS2'] <- sv[41,'POS2']
    sv[44,'POS2'] <- sv[46,'POS2']
    sv <- sv[-c(11,16:18,20,32,38:39,41,45:46,48),]
    sv <- sv[!duplicated.breakpoints(sv),]
  }
  if (gene == 'smc-5') {
    sv[16,'POS2'] <- sv[17,'POS2']
    sv <- sv[-c(4,15,17),]
  }
  if (gene == 'smc-6') {
    sv[1,'POS2'] <- sv[2,'POS2']
    sv[5,'POS2'] <- sv[7,'POS2'] 
    sv <- sv[-c(2,6:7),]
  }

  sv.properties$total_number[sv.properties$genotype == gene] <- nrow(sv)
  del.sv <- sv[sv$clust.type == 'DEL',]
  del.sv.ranges <- GRanges(seqnames = as.character(del.sv$CHR1), ranges = IRanges(start = del.sv$POS1, end = del.sv$POS2))
  td.sv <- sv[sv$clust.type == 'TD',]
  td.sv.ranges <- GRanges(seqnames = as.character(td.sv$CHR1), ranges = IRanges(start = td.sv$POS1, end = td.sv$POS2))
  
  sv.properties$svs_rep_left[sv.properties$genotype == gene] <- 0
  sv.properties$svs_rep_right[sv.properties$genotype == gene] <- 0
  if (length(td.sv.ranges)>0) {
    orientation <- list()
    for (j in 1:length(td.sv.ranges)) {
      orientation[[j]] <- c(as.character(strand(repStrand)[subjectHits(findOverlaps(query = td.sv.ranges[j], subject = repStrand))[1]]),
                            as.character(strand(repStrand)[subjectHits(findOverlaps(query = td.sv.ranges[j], subject = repStrand))[length(subjectHits(findOverlaps(query = td.sv.ranges[j], subject = repStrand)))]]))
    }
    replication_tables_td[[gene]] <- do.call('rbind',orientation)
    sv.properties$svs_rep_left[sv.properties$genotype == gene] <- sv.properties$svs_rep_left[sv.properties$genotype == gene] + 
      sum(apply(replication_tables_td[[gene]], 1, function(x) sum(x=='-')>0))
    sv.properties$svs_rep_right[sv.properties$genotype == gene] <- sv.properties$svs_rep_right[sv.properties$genotype == gene] + 
      sum(apply(replication_tables_td[[gene]], 1, function(x) sum(x=='+')>0))
  }
  
  if (length(del.sv.ranges)>0) {
    orientation <- list()
    for (j in 1:length(del.sv.ranges)) {
      orientation[[j]] <- c(as.character(strand(repStrand)[subjectHits(findOverlaps(query = del.sv.ranges[j], subject = repStrand))[1]]),
                            as.character(strand(repStrand)[subjectHits(findOverlaps(query = del.sv.ranges[j], subject = repStrand))[length(subjectHits(findOverlaps(query = del.sv.ranges[j], subject = repStrand)))]]))
    }
    replication_tables_del[[gene]] <- do.call('rbind',orientation)
    sv.properties$svs_rep_left[sv.properties$genotype == gene] <- sum(apply(replication_tables_del[[gene]], 1, function(x) sum(x=='-')>0)) + 
      sv.properties$svs_rep_left[sv.properties$genotype == gene]
    sv.properties$svs_rep_right[sv.properties$genotype == gene] <- sum(apply(replication_tables_del[[gene]], 1, function(x) sum(x=='+')>0)) + 
      sv.properties$svs_rep_right[sv.properties$genotype == gene]
  }
  
  sv.properties$no.samples[sv.properties$genotype == gene] <- sum(CD2Mutant %in% paste0(gene, endings))
  sv.properties$deletions[sv.properties$genotype == gene] <- nrow(del.sv)
  sv.properties$tandem_duplications[sv.properties$genotype == gene] <- nrow(td.sv)
  
  all.bp.ranges <- GRanges(seqnames = c(as.character(sv$CHR1),as.character(sv$CHR2)), 
                           ranges = IRanges(start = c(sv$POS1-30,sv$POS2-30), 
                                            end = c(sv$POS1+30,sv$POS2+30)))
  hits <- findOverlaps(query = all.bp.ranges, subject = repeat.table.ranges)
  sv.properties$variants_in_repeats[sv.properties$genotype == gene] <- length(unique(queryHits(hits)))
  
  hits <- findOverlaps(query = all.bp.ranges, subject = c(gc_minus,gc_plus))
  sv.properties$svs_with_Grich[sv.properties$genotype == gene] <- length(unique(queryHits(hits)))
  
  td_size[[gene]] <- width(td.sv.ranges)
  del_size[[gene]] <- width(del.sv.ranges)
  
  # indel counts - requires a list of indel VCFs indels_dedup
  #indels <- indels_dedup[names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]]
  #indels <- indels[sapply(indels,length)>0]
  #indels <- lapply(indels, granges)
  #allindels <- indels[[1]]
  #for (j in 2:length(indels))  
  #  allindels <- c(allindels,indels[[j]])
  #allindels <- unique(allindels)
  #sv.properties$indels_in_Grich[sv.properties$genotype == gene] <- 
  #  length(unique(queryHits(findOverlaps(query = allindels, subject = c(gc_minus,gc_plus))))) 
  #sv.properties$total_indels[sv.properties$genotype == gene] <- length(allindels)
  
  tel.count <- 0
  for (j in 1:nrow(sv)) {
    if (sv[j,2] < subtelomere[as.character(sv[j,1]),2] || sv[j,4] > subtelomere[as.character(sv[j,3]),3]) {
      tel.count = tel.count + 1
    }
  }
  sv.properties$variants_in_telomeres[sv.properties$genotype == gene] <- tel.count
  
  print(gene)
}

mhlist <- list()
indmhlist <- list()
for (gene in c(sv.properties$genotype, 'rev-3', 'polh-1', 'polh(lf31)-1')) {
  
  #ending = max(data$Generation[data$Genotype.new == gene])
  sv <- do.call('rbind',SVclust.new[names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]])
  sv <- sv[sv$clust.type!='some',]
  all.bp.ranges <- GRanges(seqnames = c(as.character(sv$CHR1),as.character(sv$CHR2)), 
                           ranges = IRanges(start = c(sv$POS1-30,sv$POS2-30), 
                                            end = c(sv$POS1+30,sv$POS2+30)))
  
  indels <- indels_dedup[names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]]
  indels <- indels[sapply(indels,length)>0]
  indels <- lapply(indels, granges)
  allindels <- indels[[1]]
  for (j in 2:length(indels))  
    allindels <- c(allindels,indels[[j]])
  allindels <- allindels[nchar(unlist(allindels$ALT))>5 | nchar(allindels$REF) > 5] 
  allindels <- unique(allindels)

  alist = NULL
  blist = NULL
  if (length(sv) > 0)
     for (y in 1:(length(all.bp.ranges)/2)) {
       x1 <- all.bp.ranges[y]
       x2 <- all.bp.ranges[y + length(all.bp.ranges)/2]
       seq1 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x1))),
                                           start(x1) + 1, start(x1) + 30)))
       seq2 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x2))),
                                           start(x2) + 1, start(x2) + 30)))
       i <- 30
       err <- 1
       a <- 0
       while((i>0) && err > 0) {
         if (substr(seq1,i,i) == substr(seq2,i,i)) {
           a <- a+1
           i <- i-1
         }
         else {
           err <- err - 1
         }
       }
   
       seq1 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x1))),
                                           end(x1) - 30, end(x1))))
       seq2 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x2))),
                                           end(x2) - 30, end(x2))))
       i <- 2
       err <- 1
       b <- 0
       while((i<31) && err > 0) {
         if (substr(seq1,i,i) == substr(seq2,i,i)) {
           b <- b+1
           i <- i+1
         }
         else {
           err <- err-1
         }
       }
       alist <- c(alist,a+b)
     }
  if (length(allindels)>0)
    for (y in 1:(length(allindels))) {
      if (nchar(allindels$REF)[y] < 5) {
        x <- allindels[y]
        inserted_seq <- DNAString(x = as.character(unlist(x$ALT)))
        prev_seq <- getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),start(x)-nchar(inserted_seq),start(x)-1)
        next_seq <- getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),end(x)+1, end(x)+nchar(inserted_seq))
        a = 0
        for (i in 1:nchar(inserted_seq)) {
          if (substr(inserted_seq,i,i) == substr(next_seq,i,i)) a = a + 1
          else break
        }
        for (i in nchar(inserted_seq):1) {
          if (substr(inserted_seq,i,i) == substr(prev_seq,i,i)) a = a + 1
          else break
        }
        
      } else {
        x <- allindels[y]
        deleted_seq <- DNAString(x = as.character(x$REF))
        prev_seq <- getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),start(x)-nchar(deleted_seq),start(x)-1)
        next_seq <- getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),end(x)+1, end(x)+nchar(deleted_seq))
        a = 0
        for (i in 1:nchar(deleted_seq)) {
          if (substr(deleted_seq,i,i) == substr(next_seq,i,i)) a = a + 1
          else break
        }
        for (i in nchar(deleted_seq):1) {
          if (substr(deleted_seq,i,i) == substr(prev_seq,i,i)) a = a + 1
          else break
        }
        
      }
      seq1 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),
                                          start(x) - 29, start(x))))
      seq2 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),
                                          end(x) - 29, end(x))))
      i <- 30
      err <- 1
      a <- 0
      while((i>0) && err > 0) {
        if (substr(seq1,i,i) == substr(seq2,i,i)) {
          a <- a+1
          i <- i-1
        }
        else {
          err <- err - 1
        }
      }

      seq1 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),
                                          start(x), start(x) + 30)))
      seq2 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),
                                          end(x), end(x) + 30)))
      i <- 2
      err <- 1
      b <- 0
      while((i<31) && err > 0) {
        if (substr(seq1,i,i) == substr(seq2,i,i)) {
          b <- b+1
          i <- i+1
        }
        else {
          err <- err-1
        }
      }
      blist <- c(blist,a+b)
      
    }

  mhlist[[gene]] <- alist
  indmhlist[[gene]] <- blist
  print(gene)
}  

#############################visualize######################################################

par(mar = c(6,4,2,2))  
boxplot(mhlist, las = 2)
plot(sapply(mhlist, function(x) sum(x>0)) / c(sv.properties$total_number,5,5), xaxt = 'n', 
     ylab = 'Proportion of SVs', xlab ='', pch = 16)
axis(side = 1, at = c(1:17), labels = names(mhlist), las = 2, font = 3, lty = 0)


pdf('MH_at_SV_bps.pdf', 7,5)
vioplot(mhlist, las = 2, frame = NA, xaxt = 'n', ylab = 'MH length, bp', yaxt = 'n', frame.plot = F, plotCentre='line', main = 'Microhomology length at SV breakpoints')
axis(side = 2, at = c(0:13), las = 2)
axis(side = 1, at = c(1:17), labels = names(mhlist), las = 2, font = 3, lty = 0)
for (j in 1:17)
  points(x = jitter(rep(j,length(mhlist[[j]])),amount = 0.3), y = mhlist[[j]], cex = 0.5, pch = 21, bg = 'gray86', col = 'black')
dev.off()


mhlist[['polh-1']] <- NA
mhlist <- mhlist[names(indmhlist)]

pdf('SV_properties_MH.pdf', 7, 5)
par(mar = c(6,4,2,2))  
plot(NA, NA, bty = 'n', xaxt = 'n', ylab = 'Proportion of SVs', xlab = '', ylim = c(0,1), xlim = c(1,19), las = 2, main = 'Proportion of SVs and indels with MH')
for (j in seq(2,19,2))
  polygon(x = c(j-0.5,j-0.5,j+0.5,j+0.5,j-0.5),
          y = c(0,1,1,0,0),
          col = 'grey86', border = NA)

sv.means <- sapply(mhlist, function(x) sum(x>0)) / sapply(mhlist, length)
sv.sds <- sqrt(sv.means * (1 - sv.means) / sapply(mhlist, length))
lows <- sv.means - qnorm(0.975) * sv.sds
lows[lows<0] <- 0
points(x = c(1:18)-0.2, y = sv.means,pch = 16,col = 'brown')
arrows(x0 = c(1:18) - 0.2, y0 = lows, y1 = sv.means + qnorm(0.975)*sv.sds,
       lwd = 0.5, col = 'brown', length=0, angle=90)

sv.means <- sapply(indmhlist, function(x) sum(x>0)) / sapply(indmhlist, length)
sv.sds <- sqrt(sv.means * (1 - sv.means) / sapply(indmhlist, length))
lows <- sv.means - qnorm(0.975) * sv.sds
lows[lows<0] <- 0
points(x = c(1:18)+0.2, y = sv.means,pch = 16,col = 'blue')
arrows(x0 = c(1:18)+0.2, y0 = lows, y1 = sv.means + qnorm(0.975)*sv.sds,
       lwd = 0.5, col = 'blue', length=0, angle=90)

axis(side = 1, at = c(1:18), las = 2, labels = names(mhlist), lty = 0, tick = F, tck = 0.01, font = 3)
legend('topleft', bty = 'n', legend = c('SVs with MH', 'indels with MH'), 
       col = c('brown', 'blue'), pch = 16)
dev.off()


pdf('MH_at_indel_bps.pdf', 7,5)
vioplot(indmhlist, las = 2, frame = NA, ylab = 'MH length, bp', xaxt = 'n', yaxt = 'n', frame.plot = F, plotCentre='line', 
        main = 'Microhomology length at indel breakpoints')
axis(side = 2, at = c(0,5,10,15,20,25), las = 2)
axis(side = 1, at = c(1:18), labels = names(indmhlist), las = 2, font = 3, lty = 0)
for (j in 1:18)
  points(x = jitter(rep(j,length(indmhlist[[j]])),amount = 0.3), y = indmhlist[[j]], cex = 0.5, pch = 21, bg = 'gray86', col = 'black')
dev.off()

plot(sapply(indmhlist, function(x) sum(x>0)) / sapply(indmhlist, length), xaxt = 'n', 
     ylab = 'Proportion of indels with MH', xlab ='')
axis(side = 1, at = c(1:18), labels = names(indmhlist), las = 2, font = 3, lty = 0)

# random expectation of MH
(sapply(0:28, function(x) (x+1) * (1/4)**x * (3/4)**2)) -> prblty


#############################other plots######################################################

plot(sv.properties$total_number, bty = 'n', xaxt = 'n', ylab = 'Number of SVs', xlab = '', pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)

plot(sv.properties$total_number / sv.properties$no.samples, bty = 'n', xaxt = 'n', ylab = 'Number of SVs per sample', xlab = '', pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)

plot(sv.properties$deletions/sv.properties$total_number, bty = 'n', xaxt = 'n', ylab = 'Proportion of deletions', xlab = '', pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)

plot(sv.properties$tandem_duplications/sv.properties$total_number, bty = 'n', xaxt = 'n', ylab = 'Proportion of TDs', xlab = '', pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)


allmean <- sum(width(c(gc_minus,gc_plus))) / sum(chr_sizes[-5])
allsd <- sqrt(sum(width(c(gc_minus,gc_plus))) / sum(chr_sizes[-5]) * (1 - sum(width(c(gc_minus,gc_plus))) / sum(chr_sizes[-5])))

pdf('Grich_SV_indel.pdf', 7, 5)
plot(NA, NA, bty = 'n', xaxt = 'n', ylab = 'Proportion of variants', xlab = '', ylim = c(0,1), xlim = c(1,15), las = 2)
for (j in seq(2,15,2))
  polygon(x = c(j-0.5,j-0.5,j+0.5,j+0.5,j-0.5),
          y = c(0,1,1,0,0),
          col = 'grey86', border = NA)
sv.means <- sv.properties$svs_with_Grich/(sv.properties$total_number)
sv.sds <- sqrt(sv.means * (1 - sv.means) / sv.properties$total_number)
lows <- sv.means - qnorm(0.975) * sv.sds
lows[lows<0] <- 0
abline(h = allmean, lty = 2)
points(sv.means,pch = 16,col = 'black')
arrows(x0 = 1:15, y0 = lows, y1 = sv.means + qnorm(0.975)*sv.sds,
       lwd = 0.5, col = 'black', length=0, angle=90)

sv.means <- sv.properties$indels_in_Grich/(sv.properties$total_indels)
sv.sds <- sqrt(sv.means * (1 - sv.means) / sv.properties$total_indels)
lows <- sv.means - qnorm(0.975) * sv.sds
lows[lows<0] <- 0
points(x = c(1:15)+0.2, y = sv.means,pch = 16,col = 'brown')
arrows(x0 = c(1:15) + 0.2, y0 = lows, y1 = sv.means + qnorm(0.975)*sv.sds,
       lwd = 0.5, col = 'brown', length=0, angle=90)
legend('topright', bty = 'n', legend = c('SVs in G-rich sequences', 'indels in G-rich sequences', 'expected proportion'), 
       col = c('black', 'brown', 'black'), pch = c(16,16,16))
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)
dev.off()


pdf('SV_properties_telomeres.pdf', 7, 5)
plot(NA, NA, bty = 'n', xaxt = 'n', ylab = 'Proportion of SVs', xlab = '', ylim = c(0,1), xlim = c(1,15), las = 2)
for (j in seq(2,15,2))
  polygon(x = c(j-0.5,j-0.5,j+0.5,j+0.5,j-0.5),
          y = c(0,1,1,0,0),
          col = 'grey86', border = NA)

sv.means <- sv.properties$variants_in_telomeres/(sv.properties$total_number)
sv.sds <- sqrt(sv.means * (1 - sv.means) / sv.properties$total_number)
lows <- sv.means - qnorm(0.975) * sv.sds
lows[lows<0] <- 0
points(x = c(1:15)-0.1, y = sv.means,pch = 16,col = 'brown')
arrows(x0 = c(1:15) - 0.1, y0 = lows, y1 = sv.means + qnorm(0.975)*sv.sds,
       lwd = 0.5, col = 'brown', length=0, angle=90)
alltel <- (sum(subtelomere[,4] - subtelomere[,3]) + sum(subtelomere[,2] - subtelomere[,1])) / sum(chr_sizes[-5])
abline(h = alltel / sum(chr_sizes[-5]), col = 'brown', lty = 2)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)
dev.off()


pdf('SV_properties_repeats.pdf', 7, 5)
plot(NA, NA, bty = 'n', xaxt = 'n', ylab = 'Proportion of SVs', xlab = '', ylim = c(0,1), xlim = c(1,15), las = 2)
for (j in seq(2,15,2))
  polygon(x = c(j-0.5,j-0.5,j+0.5,j+0.5,j-0.5),
          y = c(0,1,1,0,0),
          col = 'grey86', border = NA)

sv.means <- sv.properties$variants_in_repeats/(sv.properties$total_number)
sv.sds <- sqrt(sv.means * (1 - sv.means) / sv.properties$total_number)
lows <- sv.means - qnorm(0.975) * sv.sds
lows[lows<0] <- 0
points(x = c(1:15)+0.1, y = sv.means,pch = 16,col = 'blue')
arrows(x0 = c(1:15)+0.1, y0 = lows, y1 = sv.means + qnorm(0.975)*sv.sds,
       lwd = 0.5, col = 'blue', length=0, angle=90)
abline(h = sum(width(repeat.table.ranges)) / sum(chr_sizes[-5]), col = 'blue', lty = 2)
legend('topright', bty = 'n', legend = c('in telomeric regions', 'in repetitive regions', 'expected proportion'), 
       col = c( 'brown', 'blue', 'black'), pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)
dev.off()


############################################################

# plots of TD and DEL size distribution

# tandem duplications
log_td_size <- lapply(td_size, log10)
names(log_td_size) <- names(td_size)
pdf('new_sizes_of_TDs_Feb2020.pdf',8,6)
par(mar = c(9,6,2,2))
vioplot(x = log_td_size, col = brewer.pal(7,'Pastel1')[1], frame = NA, xaxt = 'n', yaxt = 'n', frame.plot = F, ylim = c(2,7), plotCentre='line')
for (j in 1:length(td_size))
  points(x = jitter(rep(j, length(log_td_size[[j]])), amount = .3), y = log_td_size[[j]], pch = 16, col= 'grey')
axis(side = 1, las = 2, labels = print_names[sapply(names(td_size), function(x) grep(x,names(print_names))[1])], font = 3, at = 1:length(td_size), lty = 0)
axis(side = 2, las = 2, at = c(2:7), labels = c('100 bps', '1 kbps', '10 kbps', '100 kbps', '1 Mbps', '10 Mbps'))
dev.off()

# deletions
del_size <- del_size[sapply(del_size, length)>0]
log_td_size <- lapply(del_size, log10)
names(log_td_size) <- names(del_size)
pdf('new_sizes_of_DELs_Feb2020.pdf',8,6)
par(mar = c(9,6,2,2))
vioplot(x = log_td_size, col = brewer.pal(7,'Pastel1')[3], frame = NA, xaxt = 'n', yaxt = 'n', frame.plot = F, ylim = c(2,6), plotCentre='line')
for (j in 1:length(del_size)) {
  if (length(log_td_size) > 1)
    points(x = jitter(rep(j, length(log_td_size[[j]])), amount = .3), y = log_td_size[[j]], pch = 16, col= 'grey')
  else points(x = j, y = log_td_size[[j]], pch = 16, col = 'grey')
}
axis(side = 1, las = 2, labels = print_names[sapply(names(del_size), function(x) grep(x,names(print_names))[1])], font = 3, at = 1:length(del_size), lty = 0)
axis(side = 2, las = 2, at = c(2:6), labels = c('100 bps', '1 kbps', '10 kbps', '100 kbps', '1 Mbps'))
dev.off()

######## session information ########
# R version 3.6.3 (2020-02-29)
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] VariantAnnotation_1.30.1          Rsamtools_2.0.3                   SummarizedExperiment_1.14.1      
# [4] DelayedArray_0.10.0               BiocParallel_1.18.1               matrixStats_0.55.0               
# [7] Biobase_2.44.0                    BSgenome.Celegans.UCSC.ce11_1.4.2 BSgenome_1.52.0                  
# [10] Biostrings_2.52.0                 XVector_0.24.0                    rtracklayer_1.44.4               
# [13] GenomicRanges_1.36.1              GenomeInfoDb_1.20.0               IRanges_2.18.3                   
# [16] S4Vectors_0.22.1                  BiocGenerics_0.30.0               reshape2_1.4.3                   
# [19] ggplot2_3.2.1                    

# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.2               lattice_0.20-38          prettyunits_1.0.2        assertthat_0.2.1         digest_0.6.22           
# [6] R6_2.4.0                 plyr_1.8.4               RSQLite_2.1.2            evaluate_0.14            httr_1.4.1              
# [11] pillar_1.4.2             progress_1.2.2           zlibbioc_1.30.0          rlang_0.4.5              GenomicFeatures_1.36.4  
# [16] lazyeval_0.2.2           rstudioapi_0.10          blob_1.2.0               Matrix_1.2-18            rmarkdown_1.16          
# [21] stringr_1.4.0            biomaRt_2.40.5           RCurl_1.95-4.12          bit_1.1-14               munsell_0.5.0           
# [26] compiler_3.6.3           xfun_0.10                pkgconfig_2.0.3          htmltools_0.4.0          tidyselect_0.2.5        
# [31] tibble_2.1.3             GenomeInfoDbData_1.2.1   XML_3.98-1.20            crayon_1.3.4             dplyr_0.8.3             
# [36] withr_2.1.2              GenomicAlignments_1.20.1 bitops_1.0-6             grid_3.6.3               gtable_0.3.0            
# [41] DBI_1.0.0                magrittr_1.5             scales_1.0.0             stringi_1.4.3            vctrs_0.2.4             
# [46] tools_3.6.3              bit64_0.9-7              glue_1.3.1               purrr_0.3.3              hms_0.5.3               
# [51] yaml_2.2.0               AnnotationDbi_1.46.1     colorspace_1.4-1         memoise_1.1.0            knitr_1.25 

