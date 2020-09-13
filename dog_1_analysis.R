#######################################################################
# a script to generate a distribution of variants in dog-1 mutations  #
# N. Volkova, EMBL-EBI, 2019                                          #
# for Meier et al. publication on DNA repair signatures in C. elegans #
#######################################################################

library(rtracklayer)
library(VariantAnnotation)

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

# load mutations
PATHTOSUBSVCF='/path/to/filtered/substitution/vcfs'
PATHTOINDELSVCF='/path/to/filtered/indel/vcfs'
PATHTOSV='/path/to/filtered/SV/tables'
indels_dedup <- sapply(data$Sample[data$Genotype.new == 'dog-1'], function(x) readVcf(paste0(PATHTOSUBSVCF,'/',x,'.vcf')))
vcfs_dedup <- sapply(data$Sample[data$Genotype.new == 'dog-1'], function(x) readVcf(paste0(PATHTOINDELSVCF,'/',x,'.vcf')))
SVclust.new <- sapply(data$Sample[data$Genotype.new == 'dog-1'], function(x) readVcf(paste0(PATHTOSV,'/',x,'.vcf')))

#####################################################################################

# load GC tracks from Marsico et al. 2019
gc_plus <- import('~/Downloads/Celegans_all_w15_th-1_plus.hits.max.PDS.w50.35.bed')
gc_minus <- import('~/Downloads/Celegans_all_w15_th-1_minus.hits.max.PDS.w50.35.bed')
seqlevels(gc_plus) <- sapply(seqlevels(gc_plus), function(x) unlist(strsplit(x,split = '[_]'))[2])
seqlevels(gc_minus) <- sapply(seqlevels(gc_minus), function(x) unlist(strsplit(x,split = '[_]'))[2])

#####################################################################################

intnames <- names(CD2Mutant)[CD2Mutant == 'dog-1:20']
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
