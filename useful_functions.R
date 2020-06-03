################################################################################
##### Supplementary functions                                              #####
##### January 2016 - ... , EMBL-EBI, N.Volkova (nvolkova@ebi.ac.uk)        #####
################################################################################

round_dataframe <- function(df, n = 2) {
  
  inds <- which(sapply(df,is.numeric))
  
  for (j in inds) df[,j] <- round(df[,j],n)
  
  return(df)
  
}

lenght <- function(x) {
  message("Learn to spell it properly!")
  return(length(x))
}

top <- function(x,n) {
  return(sort(x,decreasing=T)[1:n])
}
which.top <- function(x,n) {
  if (is.null(names(x))) return(order(x,decreasing=T)[1:n])
  return(names(x)[order(x,decreasing=T)][1:n])
}

bottom <- function(x,n) {
  return(sort(x,decreasing=F)[1:n])
}
which.bottom <- function(x,n) {
  if (is.null(names(x))) return(order(x,decreasing=F)[1:n])
  return(names(x)[order(x,decreasing=F)][1:n])
}

cosine <- function(x,y) {
    return(sum(x * y) / sqrt(sum(x**2)) / sqrt(sum(y**2)))
}

divergence <- function (a,b) {
  return (sum(a * log ( (a+.Machine$double.eps)/(b + .Machine$double.eps)) - a + b))
}

KL <- function(p,q) {
  return(sum(p * log2((p+0.001)/(q+0.001))))
}

JSdistance <- function(p,q) {
  M = 0.5*(p+q)
  return (sqrt(0.5*KL(p,M) + 0.5*KL(q,M)))
}

restore <- function(x) {
  t1 <- unlist(strsplit(x, split = '[:]'))
  t2 <- c(t1[1],unlist(strsplit(t1[2], split = '[_]')))
  t3 <- c(t2[1:2],unlist(strsplit(t2[3], split = '[/]')))
  return(t3)
}

read_ce_vcf <- function(file, genome = "WBcel235") {
  out <- tryCatch(
    {
      readVcf(file, genome=genome) 
    },
    error=function(cond) {
      message(paste("\nBad file in: ", file))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("\nFile caused a warning: ", file))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    }
  )    
  return(out)
}

combination_check <- function(x,y) {
  objective <- function(l) {
    return(-as.numeric(cosine(y,x %*% t(t(l)))))
  }
  return(optim(par = rep(1,ncol(x)), objective, lower=rep(0,ncol(x))))
}

divergence <- function (a,b) {
  return (a * log ( (a+.Machine$double.eps)/(b + .Machine$double.eps)) - a + b)
}

# supplementary functions
poisl <- function(X, theta, Y){
  Yp <- X %*% theta
  Y %*% log(Yp) - sum(Yp) - sum(log(factorial(Y)))
}
poisdl <- function(X, theta, Y){
  lambda <- as.numeric(X %*% theta)
  t(Y / lambda - 1) %*% X
}
poisI <- function(X, theta, Y){
  lambda <- as.numeric(X %*% theta)
  t(Y/lambda^2 * X) %*% X
}
poisL <- function(X, theta, Y){
  Yp <- X %*% theta
  prod(exp(-Yp)) * prod(Yp**Y) / prod(factorial(Y))
}

RevCom <- function(x) {
  return(as.character(reverseComplement(DNAString(x))))
}

plot_all_changes <- function(plot_fc,S_est,beta_est,n = F) {
  to.show <- data.frame(t(S_est))
  for (j in 1:nrow(beta_est)) {
    to.show <- cbind(to.show, S_est[nrow(S_est),] * exp(beta_est[j,]))
  }
  colnames(to.show) <- c(paste0('Sig',1:nrow(S_est)), paste0('Sig*beta',1:nrow(beta_est)))
  plot_fc(to.show, norm = n)
}

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)), col = "black", sep = 0.5, ...)
{
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
}

getTrinucleotideSubs <- function(vcf, ref_genome="BSgenome.Celegans.UCSC.ce11") {
  seqlevels(vcf) <- seqnames(get(ref_genome))
  if (ref_genome==worm_ref_genome)
    seqlevels(vcf) <- c("chrM","chrIV","chrIII","chrX","chrI","chrV","chrII")
  vcf <- vcf[as.character(seqnames(vcf)) != "chrM"]
  ranges = resize(vcf, 3, fix = "center")
  tnc = as.character(getSeq(get(ref_genome), seqnames(vcf), 
                            start(vcf) - 1, end(vcf) + 1))
  s <- paste0(substr(tnc,1,1),"[",as.character(ref(vcf)), ">",
              as.character(unlist(alt(vcf))),"]", substr(tnc,3,3))
  n <- c("A","C","G","T")
  f <- paste0(rep(n, each=4), "[", rep(n, each=96/2), ">", c(rep(c("C","G","T"), each=48/3),rep(c("A","G","T"), each=48/3),rep(c("A","C","T"), each=48/3), rep(c("A","C","G"), each=48/3)), "]", n)
  s <- factor(s, levels=f)
  return(s)
}

getTrinucleotideSubs_human_table <- function(subs, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19', chr = 0) {
  g = get(ref_genome)
  if (is.na(subs)) return(NA)
  if (chr) subs <- subs[subs$Chrom != "chrM",] else subs <- subs[subs$Chrom != "MT",]
  if (nrow(subs)==0) return(NA)
  if (chr)
    tnc = as.character(getSeq(g, subs$Chrom, subs$Pos - 1, subs$Pos + 1))
  else
    tnc = as.character(getSeq(g, paste('chr',subs$Chrom,sep=''), 
                            subs$Pos - 1, subs$Pos + 1))
  s <- paste0(substr(tnc,1,1),"[",as.character(subs$Ref), ">",
              as.character(subs$Alt),"]", substr(tnc,3,3))
  n <- c("A","C","G","T")
  f <- paste0(rep(n, each=4), "[", 
              rep(n, each=96/2), ">", 
              c(rep(c("C","G","T"), each=16),
                rep(c("A","G","T"), each=16),
                rep(c("A","C","T"), each=16), 
                rep(c("A","C","G"), each=16)),
              "]", n)
  s <- factor(s, levels=f)
  return(s)
}

getPentanucleotideSubs_human_table <- function(subs, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19') {
  g = get(ref_genome)
  if (is.na(subs)) return(NA)
  subs <- subs[subs$Chrom != "MT",]
  if (nrow(subs)==0) return(NA)
  tnc = as.character(getSeq(g, paste('chr',subs$Chrom,sep=''), 
                            subs$Pos - 2, subs$Pos + 2))
  s <- paste0(substr(tnc,1,2),"[",as.character(subs$Ref), ">",
              as.character(subs$Alt),"]", substr(tnc,4,5))
  n <- c("A","C","G","T")
  f <- paste0(rep(n, each=64), 
              rep(n, each=16),
              "[", 
              rep(n, each=1536/2),
              ">", 
              c(rep(c("C","G","T"), each=256),
                rep(c("A","G","T"), each=256),
                rep(c("A","C","T"), each=256), 
                rep(c("A","C","G"), each=256)),
              "]", rep(n, each = 4), n)
  s <- factor(s, levels=f)
  return(s)
}

tncToPyrimidine <- function(nucl) {
  ind <- sort(c(grep('[G',nucl,fixed = T),grep('[A',nucl,fixed = T)))
  newnucl <- as.character(reverseComplement(DNAStringSet(paste(substr(nucl[ind],1,1),
                                                               substr(nucl[ind],5,5),
                                                               substr(nucl[ind],3,3),
                                                               substr(nucl[ind],7,7),sep=''))))
  nucl[ind] <- paste(substr(newnucl,1,1),'[',
                     substr(newnucl,2,2),'>',
                     substr(newnucl,3,3),']',
                     substr(newnucl,4,4),sep='')
  return(nucl)
}

RevComMutType <- function(nucl) {
  tmp <- unlist(strsplit(nucl, split = ''))
  return(paste0(RevCom(tmp[7]),'[',RevCom(tmp[3]),'>',RevCom(tmp[5]),']',RevCom(tmp[1])))
}

pncToPyrimidine <- function(nucl) {
  ind <- sort(c(grep('[G',nucl,fixed = T),grep('[A',nucl,fixed = T)))
  newnucl <- as.character(reverseComplement(DNAStringSet(paste(substr(nucl[ind],1,2),
                                                               substr(nucl[ind],6,6),
                                                               substr(nucl[ind],4,4),
                                                               substr(nucl[ind],8,9),sep=''))))
  nucl[ind] <- paste(substr(newnucl,1,2),'[',
                     substr(newnucl,3,3),'>',
                     substr(newnucl,4,4),']',
                     substr(newnucl,5,6),sep='')
  return(nucl)
}

get_DNV_type <- function(vcf, ref_genome="BSgenome.Celegans.UCSC.ce11", dnv.types) {
  s <- paste0(paste(as.character(ref(vcf)),collapse=''), ">", paste(as.character(unlist(alt(vcf))), collapse=''))
  if (!(s %in% dnv.types)) {
    s <- paste0(as.character(reverseComplement(DNAString(paste(as.character(ref(vcf)),collapse='')))), ">",
                as.character(reverseComplement(DNAString(paste(as.character(unlist(alt(vcf))), collapse='')))))
  }
  return(s)
}

isMNV <- function(vcf) {
  d <- diff(start(vcf)) == 1 & abs(diff(geno(vcf)$PM[,"TUMOUR"] )) <= 0.05
  w <- c(FALSE, d) | c(d, FALSE)
  return(w)
}

is.in.repeat <- function(vcf, ref_genome="BSgenome.Celegans.UCSC.ce11") {
  seqlevels(vcf) <- c("chrM","chrIV","chrIII","chrX","chrI","chrV","chrII")
  vcf <- vcf[as.character(seqnames(vcf)) != "chrM"]
  #ranges = resize(vcf, 3, fix = "center")
  in.repeat <- sapply(1:length(vcf), function(j) {
    tmp <- vcf[j]
    if (info(tmp)$PC=="DI") return(FALSE)
    tnc = as.character(getSeq(get(ref_genome), seqnames(tmp), 
                              start(tmp) - info(vcf)$LEN[j]*5,
                              end(tmp) + info(vcf)$LEN[j]*5))
    if (info(tmp)$PC=="I")
      return(grepl(paste0('(',substr(as.character(unlist(alt(tmp))), 2, nchar(as.character(unlist(alt(tmp))))),'){3,}'), tnc))
    if (info(tmp)$PC=="D")
      return(grepl(paste0('(',substr(as.character(ref(tmp)), 2, nchar(as.character(ref(tmp)))),'){3,}'), tnc))
    return(FALSE)
  })
  return(in.repeat)
}

min_ranges <- function(gr) {
  if (length(gr) == 1) return(gr)
  to.merge <- 1
  while (sum(to.merge)!=0) {
    gr_start <- start(gr)
    gr_end <- end(gr)
    m = 1
    to.merge <- vector('numeric',length(gr))
    for (j in 2:length(gr)) {
      if (gr_start[j] < gr_end[j-1]) {
        if (gr_end[j] <= gr_end[j-1]) to.merge[j] <- -1
        else { 
          to.merge[j] <- m
          to.merge[j-1] <- m
          m <- m+1
        }
      }
    }
    if (sum(to.merge==-1)>0) {
      gr <- gr[-which(to.merge == -1)]
      to.merge <- to.merge[-which(to.merge == -1)]
    }
    if (sum(to.merge>0)>0) {
      for (j in 1:max(to.merge)) {
        start(gr[to.merge == j]) <- min(gr_start[to.merge == j])
        end(gr[to.merge == j]) <- min(gr_end[to.merge == j])
      }
    }
    gr <- unique(gr)
  }
  
  return(gr)
}

`%over.50%` <- function(query, subject) overlapsAny(query, subject, minoverlap = 0.5 * min(width(query),width(subject)))
`%over.90%` <- function(query, subject) overlapsAny(query, subject, minoverlap = 0.9 * min(width(query),width(subject)))

duplicated.breakpoints <- function(SV) {
  
  dupl <- rep(FALSE, nrow(SV))
  
  for (j in 2:nrow(SV)) {
    
    for (k in 1:(j-1)) {
      if (SV$CHR1[j] == SV$CHR1[k] && SV$CHR2[j] == SV$CHR2[k] && SV$clust.type[j] == SV$clust.type[k]) {
        if (SV$CHR1[j] != SV$CHR2[j]) {
          iranges1 <- IRanges(start=c(as.numeric(SV$POS1)[j] - 50,as.numeric(SV$POS2)[j] - 50),
                              end=c(as.numeric(SV$POS1)[j] + 50,as.numeric(SV$POS2)[j] + 50))
          iranges2 <- IRanges(start=c(as.numeric(SV$POS1)[k] - 50,as.numeric(SV$POS2)[k] - 50),
                              end=c(as.numeric(SV$POS1)[k] + 50,as.numeric(SV$POS2)[k] + 50))
          if (iranges1[1] %over.90% iranges2[1] && iranges1[2] %over.90% iranges2[2]) {
            dupl[j] <- TRUE
            break
          }
        } else {
          iranges <- IRanges(start=as.numeric(SV$POS1)[c(j,k)],end=as.numeric(SV$POS2)[c(j,k)])
          if (iranges[1] %over.90% iranges[2]) {
            dupl[j] <- TRUE
            break
          }
        }
      }
    }

  }
  
  return(dupl)
}

filter.breakpoints <- function(SV, reference) {
  
  reference$CHR1 <- as.character(reference$CHR1)
  reference$CHR2 <- as.character(reference$CHR2)
  reference$TYPE <- as.character(reference$TYPE)
  
  if (is.na(SV) || nrow(SV)==0) return(SV)
  
  bad <- NULL
  
  SV$CHR1 <- as.character(SV$CHR1)
  SV$CHR1 <- as.character(SV$CHR1)
  SV$TYPE <- as.character(SV$TYPE)
  
  `%over.50%` <- function(query, subject) overlapsAny(query, subject, minoverlap = 0.5 * min(width(query),width(subject)))
  
  for (j in 1:nrow(SV)) {
    
    tmp <- reference[reference$CHR1 == SV$CHR1[j] & 
                       reference$CHR2 == SV$CHR2[j] & 
                       reference$TYPE == SV$TYPE[j],,drop = F]
    if (nrow(tmp) > 0) {
      big.ranges.1 <- IRanges(start=as.numeric(tmp$POS1)-200,end=as.numeric(tmp$POS1)+200)
      big.ranges.2 <- IRanges(start=as.numeric(tmp$POS2)-200,end=as.numeric(tmp$POS2)+200)
      inds <- intersect(which(big.ranges.1 %over.50% IRanges(start=as.numeric(SV$POS1[j])-200,end=as.numeric(SV$POS1[j])+200)),
                        which(big.ranges.2 %over.50% IRanges(start=as.numeric(SV$POS2[j])-200,end=as.numeric(SV$POS2[j])+200)))
      if (length(inds) > 0)
        bad <- c(bad,j)
    }
  }
  if (length(bad)>0)
    return(SV[-bad,,drop=F])
  return(SV)
}

GG_save_pdf = function(list, filename, ...) {
  #start pdf
  pdf(filename, ...)
  #loop
  count <- 1
  for (p in list) {
    print(p)
    print(count)
    count = count +1
  }
  #end pdf
  dev.off()
  invisible(NULL)
}

