library(rtracklayer)

gc_plus <- import('~/Downloads/Celegans_all_w15_th-1_plus.hits.max.PDS.w50.35.bed')
gc_minus <- import('~/Downloads/Celegans_all_w15_th-1_minus.hits.max.PDS.w50.35.bed')
seqlevels(gc_plus) <- sapply(seqlevels(gc_plus), function(x) unlist(strsplit(x,split = '[_]'))[2])
seqlevels(gc_minus) <- sapply(seqlevels(gc_minus), function(x) unlist(strsplit(x,split = '[_]'))[2])

# intersect dog-1 variants with it

lapply(indels_dedup[names(CD2Mutant)[CD2Mutant == 'brc-1:20']], granges) -> a

allindels <- c(a[[1]],a[[2]],a[[3]],a[[4]],a[[5]],a[[6]])

lapply(vcfs_dedup[names(CD2Mutant)[CD2Mutant == 'brc-1:20']], granges) -> a

allsubs <-  c(a[[1]],a[[2]],a[[3]],a[[4]],a[[5]],a[[6]])

#allindels <- allindels[info(allindels)$LEN>=50]

replication_tables_td <- list()
replication_tables_del <- list()
endings <- c(":0",":1",":5",":6",":10",":20",":15",":40")

for (gene in names(print_names)) {

  sv <- do.call('rbind',SVclust.new[names(CD2Mutant)[CD2Mutant %in% paste0(gene, endings)]])
  del.sv <- sv[sv$clust.type == 'DEL',]
  del.sv.ranges <- GRanges(seqnames = del.sv$CHR1, ranges = IRanges(start = del.sv$POS1, end = del.sv$POS2))
  td.sv <- sv[sv$clust.type == 'TD',]
  td.sv.ranges <- GRanges(seqnames = td.sv$CHR1, ranges = IRanges(start = td.sv$POS1, end = td.sv$POS2))
  
  if (length(td.sv.ranges)>0) {
    orientation <- list()
    for (j in 1:length(td.sv.ranges)) {
      orientation[[j]] <- c(as.character(strand(repStrand)[subjectHits(findOverlaps(query = td.sv.ranges[j], subject = repStrand))[1]]),
                            as.character(strand(repStrand)[subjectHits(findOverlaps(query = td.sv.ranges[j], subject = repStrand))[length(subjectHits(findOverlaps(query = td.sv.ranges[j], subject = repStrand)))]]))
    }
    replication_tables_td[[gene]] <- do.call('rbind',orientation)
  }
  
  if (length(del.sv.ranges)>0) {
    orientation <- list()
    for (j in 1:length(del.sv.ranges)) {
      orientation[[j]] <- c(as.character(strand(repStrand)[subjectHits(findOverlaps(query = del.sv.ranges[j], subject = repStrand))[1]]),
                            as.character(strand(repStrand)[subjectHits(findOverlaps(query = del.sv.ranges[j], subject = repStrand))[length(subjectHits(findOverlaps(query = del.sv.ranges[j], subject = repStrand)))]]))
    }
    replication_tables_del[[gene]] <- do.call('rbind',orientation)
  }

  print(gene)
}

pv <- list()
for (gene in names(replication_tables_td)) {
  
  tmp <- c(0,0,0)
  names(tmp) <- c('+','-','*')
  tmp[names(table(replication_tables_td[[gene]]))] <-  table(replication_tables_td[[gene]])
  
  if (sum(tmp['+'] + tmp['-'])>0) {
    if (tmp['+']>0)
      pv[[gene]] <- prop.test(tmp['+'], n = tmp['+'] + tmp['-'], p = 0.5)$p.value
    else 
      pv[[gene]] <- prop.test(tmp['-'], n = tmp['+'] + tmp['-'], p = 0.5)$p.value
  }
}
which(unlist(pv)<0.05)
which(p.adjust(unlist(pv), method = 'BH')<0.05)

pv <- list()
for (gene in names(replication_tables_del)) {
  
  tmp <- c(0,0,0)
  names(tmp) <- c('+','-','*')
  tmp[names(table(replication_tables_del[[gene]]))] <-  table(replication_tables_del[[gene]])
  
  if (sum(tmp['+'] + tmp['-'])>0) {
    if (tmp['+']>0)
      pv[[gene]] <- prop.test(tmp['+'], n = tmp['+'] + tmp['-'], p = 0.5)$p.value
    else 
      pv[[gene]] <- prop.test(tmp['-'], n = tmp['+'] + tmp['-'], p = 0.5)$p.value
  }
}
which(unlist(pv)<0.05)
which(p.adjust(unlist(pv), method = 'BH')<0.05)



#########################
hits_plus <- findOverlaps(query = c(granges(allindels),granges(allsubs),sv), subject = gc_plus)
hits_minus <- findOverlaps(query = c(granges(allindels),granges(allsubs),sv), subject = gc_minus)
length(unique(c(queryHits(hits_plus),queryHits(hits_minus))))
# 113 out of 181 total
# 109 out of 139 longer than 50 bp
hits_plus <- findOverlaps(query = sv, subject = gc_plus)
hits_minus <- findOverlaps(query = sv, subject = gc_minus)
length(unique(c(queryHits(hits_plus),queryHits(hits_minus))))
# 17 out of 21

#########################

intnames <- names(CD2Mutant)[CD2Mutant == 'dog-1:20']
#intnames <- names(CD2Mutant)[grep('brc-1,cep-1',CD2Mutant)]
clrs <- c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","brown","#8DD3C7","#FFFFB3","#BEBADA")
names(clrs) <- c('C>A','C>G','C>T','T>A','T>C','T>G','DNV','D','DI','I')

tot <- list()

for (n in intnames) {
  
  tmp_vcf <- vcfs_dedup[[n]]
  tmp_vcf_df <- as.data.frame(granges(tmp_vcf))
  if (length(tmp_vcf) > 0) {
    tmp <- isMNV(tmp_vcf)
    tmp_vcf_df$Sample <- n
    tmp_vcf_df$VAF <- geno(tmp_vcf)[['PM']][,2]
    tmp_vcf_df$Type <- NA
    if (sum(tmp)>0) {
      dnvs <- tmp_vcf_df[tmp,][seq(1,sum(tmp),2),]
      dnvs$Type <- 'DNV'
    } else {
      dnvs <- tmp_vcf_df[tmp,]
    }
    subs <- tmp_vcf_df[!tmp,]
    if (length(subs) >0) {
      subs <- subs[subs$seqnames!='MtDNA',]
      subs$Type <- paste0(subs$REF,'>',unlist(sapply(subs$ALT,as.character)))
      subs$Type[subs$Type == 'A>C'] <- 'T>G'
      subs$Type[subs$Type == 'A>G'] <- 'T>C'
      subs$Type[subs$Type == 'A>T'] <- 'T>A'
      subs$Type[subs$Type == 'G>A'] <- 'C>T'
      subs$Type[subs$Type == 'G>C'] <- 'C>G'
      subs$Type[subs$Type == 'G>T'] <- 'C>A'
    }
  } else {
    subs <- tmp_vcf_df
    dnvs <- tmp_vcf_df
  }
  
  vcf <- indels_dedup[[n]]
  indels <- as.data.frame(granges(vcf))
  if (nrow(indels) > 0) {
    indels <- indels[indels$seqnames!='MtDNA',]
    indels$Sample <- n
    indels$VAF <- (geno(vcf)$PU[,"TUMOUR"] + geno(vcf)$NU[,"TUMOUR"]) / 
      (geno(vcf)$PR[,"TUMOUR"] + geno(vcf)$NR[,"TUMOUR"]) 
    indels$Type <- ifelse(nchar(indels$REF) > nchar(unlist(sapply(indels$ALT,as.character))), yes = 'D', no = 'I')
    indels$Type[nchar(indels$REF) > 1 & nchar(unlist(sapply(indels$ALT,as.character))) > 1] <- 'DI'
  }
  
  if (nrow(subs) > 0) {
    if (nrow(indels) > 0) {
      if (nrow(dnvs) >0) tot[[n]] <- rbind(subs,indels,dnvs)
      else tot[[n]] <- rbind(subs,indels)
    } else {
      if (nrow(dnvs) >0) tot[[n]] <- rbind(subs,dnvs)
      else tot[[n]] <- subs
    }
  } else {
    if (nrow(indels) > 0) {
      if (nrow(dnvs) >0) tot[[n]] <- rbind(indels,dnvs)
      else tot[[n]] <- indels
    } else {
      if (nrow(dnvs) >0) tot[[n]] <- rbind(dnvs)
      else tot[[n]] <- NULL
    }
  }
  
  if (length(tot[[n]])>0) {
    tot[[n]]$Mode <- NA
    tot[[n]]$experiment <- CD2Mutant[n]
    if (nrow(subs)>0) 
      tot[[n]]$Mode[1:nrow(subs)] <- 'A'
    if (nrow(indels) >0)
      tot[[n]]$Mode[(nrow(subs)+1):(nrow(subs) + nrow(indels))] <- 'B'
    if (nrow(dnvs)>0)
      tot[[n]]$Mode[(nrow(subs)+nrow(indels)+1):nrow(tot[[n]])] <- 'C'
    rownames(tot[[n]]) <- NULL
    tot[[n]] <- tot[[n]][order(tot[[n]]$seqnames),]
    for (ch in levels(tot[[n]]$seqnames)) {
      tot[[n]][tot[[n]]$seqnames == ch,] <- tot[[n]][tot[[n]]$seqnames == ch,][order(tot[[n]]$start[tot[[n]]$seqnames == ch]),]
    }
    if (max(table(tot[[n]]$seqnames)) < 3) {
      tot[[n]]$clust <- 1
    } else {
      k <- isClustered(tot[[n]])
      tot[[n]]$clust <- as.character(k)
      tot[[n]]$clust[tot[[n]]$clust == 'Kataegis'] <- 5
      tot[[n]]$clust[tot[[n]]$clust != 5] <- 1
    }
  }
}

tot <- do.call('rbind',tot)

newtot <- tot[,c(1,2,3,4,13:14)]
sv$width <- sv$POS2 - sv$POS1
sv <- sv[,c(1,2,4,14,6)]
sv$MODE <- 'SV'
colnames(sv) <- colnames(newtot)
newtot <- rbind(newtot, sv)

chr_lens <- c(15072434,15279421,13783801,17493829,20924180,17718942)
names(chr_lens) <- c('I','II','III','IV','V','X')
df <- data.frame(name = names(chr_lens), length = chr_lens)

var.ranges <- GRanges(seqnames = newtot$seqnames, ranges = IRanges(start = newtot$start, end = newtot$end))
hits_plus <- findOverlaps(query = var.ranges, subject = gc_plus)
hits_minus <- findOverlaps(query = var.ranges, subject = gc_minus)
newtot$GC <- as.numeric(c(1:nrow(newtot)) %in% unique(c(queryHits(hits_plus),queryHits(hits_minus))))
newtot$GC[newtot$GC>0] <- 4
newtot$GC[newtot$GC==0] <- 1
newtot$GC <- factor(newtot$GC)

clrs <- c(clrs,"darkmagenta")
names(clrs)[11] <- 'SV'

q <- ggplot() + geom_bar(data = df, aes(x = name, y = length), stat = 'identity',fill = 'white') +
      scale_y_continuous(labels = c('0 Mb','5 Mb', '10 Mb','15 Mb','20 Mb')) +
      geom_jitter(data = newtot, aes(x = seqnames, y = start, shape = Mode, size = GC)) + # col = Type
      labs(title = paste0('Mutations in dog-1 F20 mutants (11 samples)'),x = 'Chromosome',y='Position') +
      scale_shape_discrete(labels = c('Substitutions','Indels','DNVs','SV')) +
      scale_size_manual(values = c(1,4), labels = c('non-G-rich','G-rich')) +
      #scale_color_manual(values=clrs[unique(tot$Type)]) +
      guides(size=guide_legend(title="Location"), shape = guide_legend(title='Class')) + 
      theme(title = element_text(size=10))

ggsave(plot = q, device = 'pdf', filename = '~/Mutation accumulation/Mutation_accumulation_location_plots_for_dog.1.pdf', width = 8, height = 7)
