load('~/SV_filtering_DELLY_better_clustering_2020.RData')

sv.properties <- data.frame(genotype = c('N2', 'atm-1', 'brc-1', 'brc-1,ced-3', 
                                         'brc-1,cep-1', 'dog-1', 'helq-1', 'him-6', 
                                         'him-6,ced-3', 'him-6,cep-1', 'mus-81', 'parp-2', 
                                         'rip-1', 'rtel-1', 'slx-1', 'smc-5', 'smc-6', 'wrn-1'),
                            total_number = rep(0,18),
                            deletions = rep(0,18),
                            tandem_duplications = rep(0,18),
                            #inversions = rep(0,13),
                            #interchromosomal_events = rep(0,13),
                            svs_with_insertion = c(0,0,3,NA,NA,2,4,1,NA,NA,5,NA,1,2,7,5,0,2),
                            svs_with_MH = c(1,1,7,NA,NA,5,6,1,NA,NA,16,NA,2,7,5,5,3,1),
                            svs_with_Grich = c(0,0,0,0,0,15,rep(0,12)),
                            indels_in_Grich = rep(0,18),
                            total_indels = rep(0,18),
                            svs_rep_left = rep(0,18),
                            svs_rep_right = rep(0,18),
                            variants_in_telomeres = rep(0,18),
                            variants_in_repeats = rep(0,18),
                            no.samples = rep(0,18))
sv.properties <- sv.properties[-c(10,12,18),]
sv.properties$genotype <- as.character(sv.properties$genotype)

# Go through all samples
endings <- c(":5",":6",":10",":20",":15",":40")

library(rtracklayer)
# GC rich regions
gc_plus <- import('~/Downloads/Celegans_all_w15_th-1_plus.hits.max.PDS.w50.35.bed')
gc_minus <- import('~/Downloads/Celegans_all_w15_th-1_minus.hits.max.PDS.w50.35.bed')
seqlevels(gc_plus) <- sapply(seqlevels(gc_plus), function(x) unlist(strsplit(x,split = '[_]'))[2])
seqlevels(gc_minus) <- sapply(seqlevels(gc_minus), function(x) unlist(strsplit(x,split = '[_]'))[2])

# Repetitive regions
load("~/Downloads/repStrand.RData")
replication_tables_td <- list()
replication_tables_del <- list()
library(BSgenome)
worm_ref_genome <- "BSgenome.Celegans.UCSC.ce11"
library(worm_ref_genome, character.only = TRUE)
repeat.table <- do.call('rbind', lapply(c('~/Downloads/CEmapI.gz',
                                          '~/Downloads/CEmapII.gz',
                                          '~/Downloads/CEmapIII.gz',
                                          '~/Downloads/CEmapIV.gz',
                                          '~/Downloads/CEmapV.gz',
                                          '~/Downloads/CEmapX.gz'),read.table,header=F))
head(repeat.table)
colnames(repeat.table) <- c('Chromosome', 'Start_on_chrom', 'End_on_chrom', 'Name_of_similar_in_Repbase', 'Start_RepBase','End_RepBase','Orientation','Identity')
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

# Subtelomeric regions
WBcel235 <- readDNAStringSet("~/Desktop/C. elegans WB235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa") # get worm reference genome
chr_sizes <- width(WBcel235)
names(chr_sizes) <- c("I","II","III","IV","MtDNA","V","X")
genome_size = sum(as.numeric(chr_sizes))
# telomeres
subtelomere <- matrix(1,nrow=6,ncol=4,dimnames=list(c("I","II","III","IV","V","X"),c("l.start","l.end","r.start","r.end")))
subtelomere[,1] <- rep(200,1)
subtelomere[,2] <- rep(16000,1)
subtelomere[,3] <- chr_sizes[-5] - 16000
subtelomere[,4] <- chr_sizes[-5] - 200

td_size <- list()
del_size <- list()

# Can try but rather do gene-by-gene
for (gene in sv.properties$genotype) {
  
  sv <- do.call('rbind',SVclust.new[names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]])
  no.sv <- do.call('rbind', delly.filtcounts[names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]])
  sv <- sv[sv$clust.type!='some',]
  #sv.properties$total_number[sv.properties$genotype == gene] <- sum(no.sv)
  #ending <- max(data$Generation[data$Genotype.new == gene])
  #sv <- do.call('rbind',SVclust.new[names(CD2Mutant)[CD2Mutant==paste(gene, ending, sep = ':')]])
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
    sv[1,'POS2'] <- sv[2,'POS2'] #•••••••••••
    sv[5,'POS2'] <- sv[7,'POS2'] #•••••••••••
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

  #alist = list()
  #for(y in 1:(length(all.bp.ranges)/2)) {
  #    x1 <- all.bp.ranges[y]
  #    x2 <- all.bp.ranges[y + length(all.bp.ranges)/2]
  #    seq1 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x1))),
  #                                        start(x1), end(x1))))
  #    seq2 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x2))),
  #                                        start(x2), end(x2))))
  #    res <- pairwiseAlignment(seq1, seq2, type='local', substitutionMatrix=mat, gapOpening = 1, gapExtension = 3)
  #    alist[[y]] <- nchar(res)
  #}
  #sv.properties$svs_with_MH[sv.properties$genotype == gene] <- sum(unlist(alist)>5)
  #sv.properties$svs_with_MH[sv.properties$genotype == gene] <- sum(mhlist[[gene]]>0)
  
  
  #mhlens[[gene]] <- NULL
  #for (x in names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]) {
  #  if (x %in% names(delly.vcf) & length(delly.vcf[[x]])>0)
  #    mhlens[[gene]] <- c(mhlens[[gene]],info(delly.vcf[[x]])[rownames(SVclust.new[[x]]),'HOMLEN'])
  #}
  #sv.properties$svs_with_MH[sv.properties$genotype == gene] <- sum(mhlens[[gene]]>0, na.rm = T)
  
  #cur.list <- delly.vcf[intersect(names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)], names(delly.tables))]
  #cur.list <- cur.list[!is.na(cur.list)]
  #ins.count <- 0
  #for (worm in names(cur.list)) {
  #  cur.list[[worm]] <- cur.list[[worm]][overlapsAny(cur.list[[worm]],all.bp.ranges)]
  #  if (length(cur.list[[worm]])>0)
  #    ins.count <- ins.count + length(which(info(cur.list[[worm]])$INSLEN>1))
  #}
  #sv.properties$svs_with_insertion[sv.properties$genotype == gene] <- ins.count
  
  #inslens <- NULL
  #for (x in names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]) {
  #  if (x %in% names(delly.vcf) & length(delly.vcf[[x]])>0)
  #    inslens <- c(inslens,length(which(info(delly.vcf[[x]])[rownames(SVclust.new[[x]]),'INSLEN']>0)))
  #}
  
  #sv.properties$svs_with_insertion[sv.properties$genotype == gene] <- sum(inslens)
  
  print(gene)
}

#mhlist <- list()
indmhlist <- list()
for (gene in c(sv.properties$genotype, 'rev-3', 'polh-1', 'polh(lf31)-1')) {
  
  #ending = max(data$Generation[data$Genotype.new == gene])
  sv <- do.call('rbind',SVclust.new[names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]])
  #sv <- do.call('rbind',SVclust.new[names(CD2Mutant)[CD2Mutant==paste(gene, ending, sep=':')]])
  sv <- sv[sv$clust.type!='some',]
  #all.bp.ranges <- GRanges(seqnames = c(as.character(sv$CHR1),as.character(sv$CHR2)), 
  #                         ranges = IRanges(start = c(sv$POS1-30,sv$POS2-30), 
  #                                          end = c(sv$POS1+30,sv$POS2+30)))
  
  indels <- indels_dedup[names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]]
  indels <- indels[sapply(indels,length)>0]
  indels <- lapply(indels, granges)
  allindels <- indels[[1]]
  for (j in 2:length(indels))  
    allindels <- c(allindels,indels[[j]])
  allindels <- allindels[nchar(unlist(allindels$ALT))>5 | nchar(allindels$REF) > 5] 
  allindels <- unique(allindels)

  #alist = NULL
  blist = NULL
  # if (length(sv) > 0)
  #   for (y in 1:(length(all.bp.ranges)/2)) {
  #     x1 <- all.bp.ranges[y]
  #     x2 <- all.bp.ranges[y + length(all.bp.ranges)/2]
  #     seq1 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x1))),
  #                                         start(x1) + 1, start(x1) + 30)))
  #     seq2 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x2))),
  #                                         start(x2) + 1, start(x2) + 30)))
  #     i <- 30
  #     err <- 1
  #     a <- 0
  #     while((i>0) && err > 0) {
  #       if (substr(seq1,i,i) == substr(seq2,i,i)) {
  #         a <- a+1
  #         i <- i-1
  #       }
  #       else {
  #         err <- err - 1
  #       }
  #     }
  # 
  #     seq1 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x1))),
  #                                         end(x1) - 30, end(x1))))
  #     seq2 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x2))),
  #                                         end(x2) - 30, end(x2))))
  #     i <- 2
  #     err <- 1
  #     b <- 0
  #     while((i<31) && err > 0) {
  #       if (substr(seq1,i,i) == substr(seq2,i,i)) {
  #         b <- b+1
  #         i <- i+1
  #       }
  #       else {
  #         err <- err-1
  #       }
  #     }
  #     alist <- c(alist,a+b)
  #   }
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
      #seq1 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),
      #                                    start(x) - 29, start(x))))
      #seq2 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),
      #                                    end(x) - 29, end(x))))
      #i <- 30
      #err <- 1
      #a <- 0
      #while((i>0) && err > 0) {
      #  if (substr(seq1,i,i) == substr(seq2,i,i)) {
      #    a <- a+1
      #    i <- i-1
      #  }
      #  else {
      #    err <- err - 1
      #  }
      #}

      #seq1 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),
      #                                    start(x), start(x) + 30)))
      #seq2 <- toupper(as.character(getSeq(get(worm_ref_genome), paste0('chr',as.character(seqnames(x))),
      #                                    end(x), end(x) + 30)))
      #i <- 2
      #err <- 1
      #b <- 0
      #while((i<31) && err > 0) {
      #  if (substr(seq1,i,i) == substr(seq2,i,i)) {
      #    b <- b+1
      #    i <- i+1
      #  }
      #  else {
      #    err <- err-1
      #  }
      #}
      #blist <- c(blist,a+b)
      blist <- c(blist,a)
    }

  #mhlist[[gene]] <- alist
  indmhlist[[gene]] <- blist
  print(gene)
}  
par(mar = c(6,4,2,2))  
boxplot(mhlist, las = 2)
plot(sapply(mhlist, function(x) sum(x>0)) / c(sv.properties$total_number,5,5), xaxt = 'n', 
     ylab = 'Proportion of SVs', xlab ='', pch = 16)
axis(side = 1, at = c(1:17), labels = names(mhlist), las = 2, font = 3, lty = 0)

library(vioplot)
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


#############################visualize######################################################

head(sv.properties)

plot(sv.properties$total_number, bty = 'n', xaxt = 'n', ylab = 'Number of SVs', xlab = '', pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)

plot(sv.properties$total_number / sv.properties$no.samples, bty = 'n', xaxt = 'n', ylab = 'Number of SVs per sample', xlab = '', pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)

plot(sv.properties$deletions/sv.properties$total_number, bty = 'n', xaxt = 'n', ylab = 'Proportion of deletions', xlab = '', pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)

plot(sv.properties$tandem_duplications/sv.properties$total_number, bty = 'n', xaxt = 'n', ylab = 'Proportion of TDs', xlab = '', pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)

# based on Bettina's data, could not count properly...
plot(sv.properties$svs_with_insertion/(sv.properties$total_number), bty = 'n', xaxt = 'n', 
     ylab = 'Proportion of SVs with insertion', 
     xlab = '', pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)

# what did Bettina count as microhomology???
plot(sv.properties$svs_with_MH/(sv.properties$total_number), bty = 'n', xaxt = 'n', 
     ylab = 'Proportion of SVs with MH', xlab = '',pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)

pdf('SV_properties_old.pdf', 7, 5)
par(mar = c(6,4,2,2))
plot(NA, NA, bty = 'n', xaxt = 'n', ylab = 'Proportion of SVs', xlab = '', ylim = c(0,1), xlim = c(1,16), las = 2)
for (j in seq(2,16,2))
  polygon(x = c(j-0.5,j-0.5,j+0.5,j+0.5,j-0.5),
          y = c(0,1,1,0,0),
          col = 'grey86', border = NA)
#points(x = 1:15 + 0.2, sv.properties$svs_with_insertion/(sv.properties$total_number), bty = 'n', xaxt = 'n', pch = 16, col = 'orange')
totnum <- c(3,13,15,NA,NA,30,35,19,NA,45,15,12,33,23,11,10)
sv.means <- sv.properties$svs_with_MH/(totnum)
sv.sds <- sqrt(sv.means * (1 - sv.means) / totnum)
lows <- sv.means - qnorm(0.975) * sv.sds
lows[lows<0] <- 0
arrows(x0 = 1:16, y0 = lows, y1 = sv.means + qnorm(0.975)*sv.sds,
       lwd = 0.5, col = 'purple', length=0, angle=90)
points(x = 1:16, sv.properties$svs_with_MH/totnum, pch = 16, col = 'purple')
#legend('topright', bty = 'n', legend = c('SVs with MH', 'SVs with insertions'), 
#       col = c('purple', 'orange'), pch = 16)
axis(side = 1, at = c(1:16), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)
dev.off()

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
#abline(h = allmean, lty = 2, col = 'brown')
#polygon(x = c(-0.5,-0.5,15.5,15.5,-0.5), 
#        y = c(0, 
#              allmean + 2 * allsd, 
#              allmean + 2 * allsd,
#              0, 0),
#        col = adjustcolor('brown', alpha = 0.1), border = NA)
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

plot(log10(sv.properties$svs_rep_left/sv.properties$svs_rep_right), bty = 'n', xaxt = 'n', ylab = 'Ratio of left to right replicated SVs', xlab = '', yaxt = 'n',pch = 16)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)
axis(side = 2, at = log10(c(0.5, 1, 1.5, 2,3,4,5,6)), las = 2, labels = c(0.5, 1, 1.5,2,3,4,5,6))
abline(h = 0, lty = 2)

par(mar = c(6,4,4,6))
plot(NA, NA, bty = 'n', xaxt = 'n', 
     ylab = 'Proportion of SVs with BP in repetitive regions', xlab = '',main = 'Repetitive regions',
     ylim = c(0,1), xlim = c(1,15), las = 2)
axis(side = 1, at = c(1:15), las = 2, labels = sv.properties$genotype, lty = 0, tick = F, tck = 0.01, font = 3)
for (j in 1:15)
  polygon(x = c(j - 0.4, j - 0.4, j + 0.4, j + 0.4, j + 0.4),
          y = c(0, sv.properties$total_number[j] / 100,sv.properties$total_number[j] / 100,0,0), col='brown')
points(x = 1:15, y = sv.properties$variants_in_repeats/sv.properties$total_number, pch = 16)
axis(side = 4, at = c(0,0.2,0.4,0.6,0.8,1), las = 2, lwd = 0, lwd.ticks = 2, col.text = 'brown',
     labels = c(0,20,40,60,80,100), col = 'brown')
mtext("Total number of SVs", side = 4, line = 3)
#dev.off()

tmp <- import('1018/CD0001c/CD0001c.bw')

bins <- c(seq(1,chr_lens[1],10000),
          seq(1,chr_lens[2],10000),
          seq(1,chr_lens[3],10000),
          seq(1,chr_lens[4],10000),
          seq(1,chr_lens[5],10000),
          seq(1,chr_lens[6],10000))
new_ranges <- GRanges(seqnames = rep(names(chr_lens), c(1508,1528,1379,1750,2093,1772)), ranges = IRanges(start = bins, end = bins + 100000))

scores <- NULL
for (j in 1:length(new_ranges)) {
  tmp1 <- import('1018/CD0001c/CD0001c.bw', which = new_ranges[j])
  scores <- c(scores,median(tmp1$score, na.rm = T))
  print(j)
}

scores_helq <- NULL
for (j in 1:length(new_ranges)) {
  tmp1 <- import('1018/CD0007c/CD0007c.bw', which = new_ranges[j])
  scores_helq <- c(scores_helq,median(tmp1$score, na.rm = T))
  print(j)
}

######################################################################3

load('nb5/WORM_FILTERING/DELLY_SV_after_QC_most_recent_more_filt_2020.RData')
delly.vcf <- lapply(delly.vcf, function(vcf) vcf[seqnames(vcf) %in% c("I","II","III","IV","V","X")])
delly.vcf <- lapply(delly.vcf, function(vcf) vcf[geno(vcf)[['DR']][,2] < 150 & 
                                                   geno(vcf)[['DR']][,2] > 15 &
                                                   (geno(vcf)[['DR']][,1] + geno(vcf)[['DV']][,1]) > 15 &
                                                   (geno(vcf)[['DV']][,1] / (geno(vcf)[['DR']][,1] + geno(vcf)[['DV']][,1])) > 0.05])
barplot(sapply(delly.vcf,length))
# finished filtering...
load('nb/WORM_FILTERING/DELLY_SV_after_QC_most_recent_more_filt_2020.RData')
delly.tables <- lapply(delly.vcf, function(vcf) {
  if (length(vcf)==0) return(NA)
  tmp <- data.frame(CHR1 = seqnames(granges(vcf)),
                    POS1 = start(granges(vcf)),
                    CHR2 = info(vcf)$CHR2,
                    POS2 = info(vcf)$END,
                    READS = info(vcf)$PE,
                    TYPE = info(vcf)$SVTYPE)
  rownames(tmp) <- names(vcf)
  return(tmp)
})

# upload the list with TD/DEL p-values and merge with it



############################################################

# re-do the plots of TD and DEL size distribution

# tandem duplications
log_td_size <- lapply(td_size, log10)
names(log_td_size) <- names(td_size)
pdf('new_sizes_of_TDs_Feb2020.pdf',8,6)
library(vioplot)
library(RColorBrewer)
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
library(vioplot)
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
