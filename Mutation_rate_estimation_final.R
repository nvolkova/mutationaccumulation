############################################
# Calculating mutation rates and           #
# mutation accumulation signatures         #
# in C. elegans                            #
# N. Volkova, EMBL-EBI, 2019-2020          #
# nvolkova@ebi.ac.uk                       #
############################################

# Libraries
library(greta)
library(bayesplot)
library(MASS)
library(VariantAnnotation)
library(ggplot2)
library(reshape2)
source('plotting_functions.R')
source('useful_functions.R')

##################################################################

# Upload the metadata and mutation counts

data <- read.csv("Sample_annotation_table.csv")
data$Sample <- as.character(data$Sample)
data$Genotype.new <- as.character(data$Genotype.new)
data$Code <- as.character(data$Code)
CD2Mutant <- sapply(data$Code, function(x) {
  t <- unlist(strsplit(x,split="[:]"))
  t[t=="NA"] <- ""
  return(paste(t[3],t[7],sep=":")) # genotype, experiment type, generation
}) # short notations
names(CD2Mutant) <- data$Sample

Y <- read.csv("Mutation_counts.csv")

print_names <- c('N2',
                 'agt-2',
                 'apn-1',
                 'exo-3',
                 'ndx-4',
                 'parp-1',
                 'parp-2',
                 'tdpo-1',
                 'ung-1',
                 'xpa-1',
                 'xpf-1',
                 'xpg-1',
                 'xpc-1',
                 'csb-1',
                 'polh-1 (ok3317)',
                 'polh-1 (lf31)',
                 'polk-1',
                 'rev-1 (gk147834)',
                 'rev-1 (gk924750)',
                 'rev-3',
                 'fcd-2',
                 'fnci-1',
                 'fncm-1',
                 'helq-1',
                 'dog-1',
                 'fan-1',
                 'cku-80',
                 'lig-4',
                 'polq-1',
                 'brc-1 brd-1',
                 'rip-1',
                 'smc-5',
                 'smc-6',
                 'slx-1',
                 'mus-81',
                 'smg-1',
                 'rcq-5',
                 'wrn-1',
                 'rtel-1',
                 'dna-2',
                 'bub-3 (gt2000)',
                 'bub-3 (ok3437)',
                 'san-1',
                 'atm-1',
                 'cep-1',
                 'ced-3',
                 'ced-4',
                 'brc-1 brd-1; ced-3',
                 'brc-1 brd-1; cep-1',
                 'him-6',
                 'him-6 ced-3',
                 'him-6; ced-4',	
                 'him-6; cep-1',
                 'brd-1',
                 'gen-1',
                 'lem-3',
                 'mlh-1',
                 'pms-2',
                 'pole-4',
                 'pole-4; pms-2',
                 'rfs-1',
                 'rif-1')

names(print_names) <- c(sapply(print_names, function(x) unlist(strsplit(x,split = ' '))[1])[1:7],
                        "tdp-1",
                        sapply(print_names, function(x) unlist(strsplit(x,split = ' '))[1])[9:15],
                        "polh(lf31)-1",'polk-1',
                        'rev(gk147834)-1',
                        sapply(print_names, function(x) unlist(strsplit(x,split = ' '))[1])[19:29],
                        "brc-1",
                        sapply(print_names, function(x) unlist(strsplit(x,split = ' '))[1])[31:40],
                        "bub(gt2000)-3","bub(ok3437)-3",
                        sapply(print_names, function(x) unlist(strsplit(x,split = ' '))[1])[43:47],
                        "brc-1,ced-3","brc-1,cep-1",
                        "him-6","him-6,ced-3","him-6,ced-4","him-6,cep-1",
                        sapply(print_names, function(x) unlist(strsplit(x,split = ' '))[1])[54:59],
                        "pole-4,pms-2",sapply(print_names, function(x) unlist(strsplit(x,split = ' '))[1])[61:62])

generation.function <- function(N) {
  if (is.na(N)) return(N)
  if (N==0) return(0)
  if (N==1) return(1)
  alpha = sum(sapply(1:N, function(i) 1/(2^(i-1)))) +
    0.25*sum(sapply(1:(N-1), function(i) sum(sapply(1:i, function(j) 1/(2^(j-1))))))
  return(alpha)
}

#################################################################3

# Model matrix

counts <- rowSums(Y[is.na(data$Mutagen[match(rownames(Y), data$Sample)]),]) 
generations <- sapply(data$Generation[match(names(counts), data$Sample)], generation.function)
names(generations) <- names(counts)
counts <- counts[!is.na(generations)]
generations <- generations[!is.na(generations)]
genotypes <- t(sapply(names(counts), function(x) as.numeric(unique(data$Genotype.new) == data$Genotype.new[match(x,data$Sample)])))
colnames(genotypes) <- unique(data$Genotype.new)

genotypes <- genotypes[,colSums(genotypes* generations)>0]
counts <- counts[rowSums(genotypes* generations) >0 ]
genotypes <- genotypes[rowSums(genotypes* generations)>0, ]
generations <- generations[names(counts)]

count.data <- data.frame(genotypes * generations)
ma.Y <- Y[rownames(count.data),]
count.data <- count.data[rowSums(count.data)>0,]
ma.Y <- ma.Y[rownames(count.data),]

# Estimation

#subs
counts <- rowSums(ma.Y[,1:96])
sigma <- variable(lower = 0)
rates <-  lognormal(meanlog = 0, sdlog = sigma, dim = c(1, m))
mu = (count.data %*% t(rates))
counts <- t(t(counts))
distribution(counts) = poisson(lambda = mu)
ma.model <- model(rates,sigma)
ma.draws <- mcmc(ma.model, warmup = 1000, n_samples = 1000) # do on cluster
draws_all <- do.call('rbind',ma.draws)
rates_sub <- list()
rates_sub$est<- colMeans(draws_all[,grep('rates',colnames(draws_all))])
rates_sub$low <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.025)
rates_sub$high <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.975)
rates_sub$var <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, var)

#indels
counts <- rowSums(ma.Y[,99:112])
sigma <- variable(lower = 0)
rates <-  lognormal(meanlog = 0, sdlog = sigma, dim = c(1, m))
mu = (count.data %*% t(rates))
counts <- t(t(counts))
distribution(counts) = poisson(lambda = mu)
ma.model <- model(rates,sigma)
ma.draws <- mcmc(ma.model, warmup = 1000, n_samples = 1000) # do on cluster
draws_all <- do.call('rbind',ma.draws)
rates_ind <- list()
rates_ind$est<- colMeans(draws_all[,grep('rates',colnames(draws_all))])
rates_ind$low <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.025)
rates_ind$high <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.975)
rates_ind$var <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, var)

#sv
counts <- rowSums(ma.Y[,113:119])
sigma <- variable(lower = 0)
rates <-  lognormal(meanlog = 0, sdlog = sigma, dim = c(1, m))
mu = (count.data %*% t(rates))
counts <- t(t(counts))
distribution(counts) = poisson(lambda=mu)
ma.model <- model(rates,sigma)
ma.draws <- mcmc(ma.model, warmup = 1000, n_samples = 1000) # do on cluster
draws_all <- do.call('rbind',ma.draws)
rates_sv <- list()
rates_sv$est<- colMeans(draws_all[,grep('rates',colnames(draws_all))])
rates_sv$low <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.025)
rates_sv$high <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.975)
rates_sv$var <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, var)

# full signatures
counts <- ma.Y[,1:119]
sigma <- variable(lower = 0)
rates <-  lognormal(meanlog = 0, sdlog = sigma, dim = c(ncol(ma.Y)-1, m))
mu = (count.data %*% t(rates))
size = 100
prob = size / (size + mu)
counts <- t(t(counts))
distribution(counts) = negative_binomial(prob = prob,
                                         size = size * matrix(1,nrow = nrow(count.data), ncol = ncol(ma.Y)-1))
ma.model <- model(rates,sigma)
ma.draws <- mcmc(ma.model, warmup = 2000, n_samples = 1000) # do on cluster
draws_all <- do.call('rbind',ma.draws)
beta_GH_greta <- matrix(colMeans(draws_all[,grep('rates',colnames(draws_all))]),nrow=ncol(ma.Y)-1)
beta_GH_greta_low <- matrix(apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.025),nrow=ncol(ma.Y)-1)
beta_GH_greta_high <- matrix(apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.975),nrow=ncol(ma.Y)-1)
beta_GH_greta_var <- matrix(apply(draws_all[,grep('rates',colnames(draws_all))], 2, var),nrow=ncol(ma.Y)-1)
colnames(beta_GH_greta) = colnames(beta_GH_greta_high) = colnames(beta_GH_greta_low) = colnames(beta_GH_greta_var) <- colnames(count.data)
save(beta_GH_greta,beta_GH_greta_low,beta_GH_greta_high,beta_GH_greta_var, file = 'full_beta_GH_for_MA_full.RData')

rownames(rates_sv$est) = rownames(rates_sv$low) = rownames(rates_sv$high) = rownames(rates_sv$var) <- 
  colnames(count.data)

DR <- 'agt.2'
BER <- sort(c('apn.1', 'exo.3', 'ndx.4', 'parp.1', 'parp.2','tdp.1','ung.1'))
NER <- sort(c('csb.1','xpa.1','xpc.1','xpf.1','xpg.1'))
TLS <- c('polh.1', "polh.lf31..1", "polk.1", 'polq.1', "rev.gk147834..1", 'rev.1', 'rev.3')
ICLR <- sort(c('dog.1', 'fan.1', 'fcd.2', 'fnci.1','fncm.1','helq.1'))
DSBR <- c('cku.80', 'lig.4', sort(c('brc.1', 'brd.1', 'dna.2','rif.1','rfs.1','rip.1','smc.5','smc.6',
          'mus.81','slx.1','gen.1','lem.3','smg.1','him.6','wrn.1','rcq.5','rtel.1')))
DS <- c('atm.1','ced.3','ced.4','cep.1',
             "brc.1.ced.3","brc.1.cep.1","him.6.ced.3","him.6.ced.4","him.6.cep.1")
MMR <- sort(c("pms.2","pole.4","pole.4.pms.2","mlh.1"))
SAC <- c("bub.gt2000..3","bub.ok3437..3",'san.1')
o <- match(c('N2',DR,BER,NER,DSBR,TLS,ICLR,SAC,DS,MMR), rownames(rates_sub$est))
genotype_labels <- print_names[sapply(rownames(rates_sub$est)[o], function(x) 
  grep(x,names(print_names))[1])]

final_colors <- c('#CC6677','#44AA99', '#332288')

########################################################

pdf('Mutation_rates_sub_26April2020.pdf', 8, 5)
par(mar = c(8,4,4,1))
plot(NA, NA, 
     xlim = c(0,length(genotype_labels)), 
     ylim = log10(c(0.1,110)),
     xaxt = 'n', 
     frame = F, 
     ylab = 'Het. mut-s per gen.',
     las = 2, 
     xlab = '', 
     main = 'Mutation rates', 
     yaxt = 'n')
for (j in seq(2,length(genotype_labels),2))
  polygon(x = c(j-0.3,j-0.3,j+0.7,j+0.7,j-0.3),
          y = log10(c(0.1,110,110,0.1,0.1)),
          col = 'grey86', border = NA)
start_sep = 1.7
abline(v = start_sep)
abline(v = start_sep + length(DR))
abline(v = start_sep + length(DR) + length(BER))
abline(v = start_sep + length(DR) + length(BER) + length(NER))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + length(DSBR))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + length(DSBR) + length(TLS))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + 
         length(DSBR) + length(TLS) + length(ICLR))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + 
         length(DSBR) + length(TLS) + length(ICLR) + length(SAC))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + 
         length(DSBR) + length(TLS) + length(ICLR) + length(SAC) + length(DS))
points(log10(rates_sub$est$x[o]), 
     cex = 0.5, 
     pch = 16, 
     col = final_colors[1])
axis(side = 2, 
     at = log10(c(0.1, 0.5, 1, 2, 5, 10, 50, 100)), 
     labels = c('<0.1', 0.5, 1, 2, 5, 10, 50, 100), 
     las = 2, 
     cex.axis = 0.7)

abline(h = log10(rates_sub$est['N2','x']), lty = 2, col = final_colors[1])
axis(side = 1, at = 1:nrow(rates_sub$est) + 0.2, 
     labels = genotype_labels, 
     las = 2, 
     hadj = 1, tick = F,
     cex.axis = 0.6, font = 3)
arrows(x0 = 1:nrow(rates_sub$est) , 
       y0 = log10(rates_sub$low$x[o]),
       y1 = log10(rates_sub$high$x[o]),
       col =final_colors[1], length=0, lwd=1)

points(x = c(1:62), y = log10(rates_sub$est$x[o]), pch = 16, cex = 0.5, col= final_colors[1])
pv <- list()
pv$sub <- sapply(rownames(rates_sub$est)[-1], function(z) {
  stat_mu = rates_sub$est[z,'x'] - rates_sub$est['N2','x']
  stat_sd = sqrt(rates_sub$var[z,'x'] + rates_sub$var['N2','x'])
  zscore = stat_mu / stat_sd
  return(1 - pchisq(q = zscore**2, df = 1))
})
names(which(p.adjust(pv$sub,method='BH')<0.05)) -> o2
points(x = (1:nrow(rates_sub$est))[match(o2, rownames(rates_sub$est)[o])] ,
       y = log10(rates_sub$est[o2,'x']), col = c(adjustcolor('darkred', alpha = 0.4),
                                                 'darkred')[as.numeric(rates_sub$est[o2,'x'] > 
                                                                         2 * rates_sub$est['N2','x'] |
                                                                         rates_sub$est[o2,'x'] < 
                                                                         0.5 * rates_sub$est['N2','x']) + 1], 
       pch = 16, cex = 1)
names(which(p.adjust(pv$ind,method='BH')<0.05)) -> o2

legend('topleft', legend = c('base substitutions','indels','structural variants',
                             'significantly different from N2 (FDR 5%)',
                             'significantly different from N2 (FDR 5%) and over 2-fold difference'), 
       ncol = 1, 
       col = c(final_colors,
               adjustcolor('gray', alpha = 0.4),
               'black'), 
       border = F, bty = "n", pch = 16, cex = 0.7, pt.cex = c(rep(0.5,3),rep(1,6)))
dev.off()


# indels
pdf('Mutation_rates_ind_26April2020.pdf', 8, 5)
par(mar = c(8,4,4,1))
plot(NA, NA, 
     xlim = c(0,length(genotype_labels)), 
     ylim = log10(c(0.01,110)),
     xaxt = 'n', 
     frame = F, 
     ylab = 'Het. mut-s per gen.',
     las = 2, 
     xlab = '', 
     main = 'Mutation rates', 
     yaxt = 'n')
for (j in seq(2,length(genotype_labels),2))
  polygon(x = c(j-0.3,j-0.3,j+0.7,j+0.7,j-0.3),
          y = log10(c(0.01,110,110,0.01,0.01)),
          col = 'grey86', border = NA)
start_sep = 1.7
abline(v = start_sep)
abline(v = start_sep + length(DR))
abline(v = start_sep + length(DR) + length(BER))
abline(v = start_sep + length(DR) + length(BER) + length(NER))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + length(DSBR))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + length(DSBR) + length(TLS))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + 
         length(DSBR) + length(TLS) + length(ICLR))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + 
         length(DSBR) + length(TLS) + length(ICLR) + length(SAC))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + 
         length(DSBR) + length(TLS) + length(ICLR) + length(SAC) + length(DS))

axis(side = 2, 
     at = log10(c(0.01, 0.1, 0.5, 1, 2, 5, 10, 50, 100)), 
     labels = c('<0.01',0.1, 0.5, 1, 2, 5, 10, 50, 100), 
     las = 2, 
     cex.axis = 0.7)

abline(h = log10(rates_ind$est['N2','x']), lty = 2, col = final_colors[2])

axis(side = 1, at = 1:nrow(rates_sub$est) + 0.2, 
     labels = genotype_labels, 
     las = 2, 
     hadj = 1, tick = F,
     cex.axis = 0.6, font = 3)

arrows(x0 = 1:nrow(rates_ind$est) + 0.2, 
       y0 = log10(rates_ind$low$x[o]),
       y1 = log10(rates_ind$high$x[o]),
       col =final_colors[2], length=0, lwd=1)

points(x = c(1:62) + 0.2, y = log10(rates_ind$est$x[o]), pch = 16, cex = 0.5, col = final_colors[2])

pv$ind <- sapply(rownames(rates_ind$est)[-1], function(z) {
  stat_mu = rates_ind$est[z,'x'] - rates_ind$est['N2','x']
  stat_sd = sqrt(rates_ind$var[z,'x'] + rates_ind$var['N2','x'])
  zscore = stat_mu / stat_sd
  return(1 - pchisq(q = zscore**2, df = 1))
})

names(which(p.adjust(pv$ind,method='BH')<0.05)) -> o2
points(x = (1:nrow(rates_ind$est))[match(o2, rownames(rates_ind$est)[o])] + 0.2,
       y = log10(rates_ind$est[o2,'x']), col = c(adjustcolor('darkgreen', alpha = 0.4),
                                                 'darkgreen')[as.numeric(rates_ind$est[o2,'x'] > 
                                                                           2 * rates_ind$est['N2','x'] |
                                                                           rates_ind$est[o2,'x'] < 
                                                                           0.5 * rates_ind$est['N2','x']) + 1], 
       pch = 16, cex = 1)
dev.off()

# SV
pdf('Mutation_rates_sv_26April2020.pdf', 8, 4)
par(mar = c(8,4,4,1))
plot(NA, NA, 
     xlim = c(0,length(genotype_labels)), 
     ylim = log10(c(0.01,3)),
     xaxt = 'n', 
     frame = F, 
     ylab = 'Het. mut-s per gen.',
     las = 2, 
     xlab = '', 
     main = 'Mutation rates', 
     yaxt = 'n')
for (j in seq(2,length(genotype_labels),2))
  polygon(x = c(j-0.3,j-0.3,j+0.7,j+0.7,j-0.3),
          y = log10(c(0.01,3,3,0.01,0.01)),
          col = 'grey86', border = NA)
start_sep = 1.7
abline(v = start_sep)
abline(v = start_sep + length(DR))
abline(v = start_sep + length(DR) + length(BER))
abline(v = start_sep + length(DR) + length(BER) + length(NER))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + length(DSBR))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + length(DSBR) + length(TLS))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + 
         length(DSBR) + length(TLS) + length(ICLR))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + 
         length(DSBR) + length(TLS) + length(ICLR) + length(SAC))
abline(v = start_sep + length(DR) + length(BER) + length(NER) + 
         length(DSBR) + length(TLS) + length(ICLR) + length(SAC) + length(DS))

axis(side = 2, 
     at = log10(c(0.01, 0.1, 0.5, 1, 2, 3)), 
     labels = c('<0.01',0.1, 0.5, 1, 2, 3), 
     las = 2, 
     cex.axis = 0.7)

abline(h = log10(rates_sv$est['N2','x']), lty = 2, col = final_colors[3])

axis(side = 1, at = 1:nrow(rates_sub$est) + 0.2, 
     labels = genotype_labels, 
     las = 2, 
     hadj = 1, tick = F,
     cex.axis = 0.6, font = 3)

arrows(x0 = 1:nrow(rates_sv$est) + 0.4, 
       y0 = log10(rates_sv$low$x[o]),
       y1 = log10(rates_sv$high$x[o]),
       col =final_colors[3], length=0, lwd=1)

yvalues <- rates_sv$est$x[o]; yvalues[yvalues < 0.01] <- 0.01
points(x = c(1:62) + 0.4, y = log10(yvalues), pch = 16, cex = 0.5, col = final_colors[3])

pv$sv <- sapply(rownames(rates_sv$est)[-1], function(z) {
  stat_mu = rates_sv$est[z,'x'] - rates_sv$est['N2','x']
  stat_sd = sqrt(rates_sv$var[z,'x'] + rates_sv$var['N2','x'])
  zscore = stat_mu / stat_sd
  return(1 - pchisq(q = zscore**2, df = 1))
})

names(which(p.adjust(pv$sv,method='BH')<0.05)) -> o2
points(x = (1:nrow(rates_sv$est))[match(o2, rownames(rates_sv$est)[o])] + 0.4,
       y = log10(rates_sv$est[o2,'x']), col = c(adjustcolor('darkblue', alpha = 0.4),
                                                'darkblue')[as.numeric(rates_sv$est[o2,'x'] > 
                                                                         2 * rates_sv$est['N2','x'] |
                                                                         rates_sv$est[o2,'x'] < 
                                                                         0.5 * rates_sv$est['N2','x']) + 1], 
       pch = 16, cex = 1)

dev.off()

# test if the rates are significantly different from N2

# for total rates
p.adjust(pv$sub, method = 'BH') -> rate.qvalues
p.adjust(pv$ind, method = 'BH') -> ind.qvalues
p.adjust(pv$sv, method = 'BH') -> sv.qvalues

#############################################################################

# Analyze full profiles

comparisons <- matrix(NA, nrow = nrow(beta_GH_greta), ncol = ncol(beta_GH_greta)-1, dimnames = list(colnames(Y)[1:119], colnames(beta_GH_greta)[-1]))
for (z in colnames(beta_GH_greta)[-1]) {
  stat_mu = beta_GH_greta[,z] - beta_GH_greta[,'N2']
  stat_sd = sqrt(beta_GH_greta_var[,z] + beta_GH_greta_var[,'N2'])
  zscores = stat_mu / stat_sd
  comparisons[,z] <- 1-pnorm(q=zscores, mean=0,sd=1)
}

comparisons <- matrix(p.adjust(comparisons,method='BH'), 
                      nrow = nrow(beta_GH_greta), 
                      ncol = ncol(beta_GH_greta)-1, 
                      dimnames = list(colnames(Y)[1:119], colnames(beta_GH_greta)[-1]))
unique(colnames(comparisons)[which(comparisons<0.05,arr.ind = T)[,2]])

comparisons.hl <- matrix(NA, nrow = 11, ncol = ncol(beta_GH_greta)-1, dimnames = list(c(c("C>A","C>G","C>T","T>A","T>C","T>G"), 'MNV', 'Deletions',
                                                                                        'Complex\n indels','Insertions', 'Structural\n Variants'), 
                                                                                      colnames(beta_GH_greta)[-1]))
for (z in colnames(beta_GH_greta)[-1]) {
  
  stat_mu = sum(beta_GH_greta[1:16,z]) - sum(beta_GH_greta[1:16,'N2'])
  stat_sd = sqrt(sum(beta_GH_greta_var[1:16,z] + beta_GH_greta_var[1:16,'N2']))
  zscores = stat_mu / stat_sd
  comparisons.hl[1,z] <- 1 - pchisq(q = zscores**2, df = 1)#1-pnorm(q=zscores, mean=0,sd=1)#1 - pchisq(q = zscores**2, df = 1)
  
  stat_mu = sum(beta_GH_greta[17:32,z] - beta_GH_greta[17:32,'N2'])
  stat_sd = sqrt(sum(beta_GH_greta_var[17:32,z] + beta_GH_greta_var[17:32,'N2']))
  zscores = stat_mu / stat_sd
  comparisons.hl[2,z] <- 1 - pchisq(q = zscores**2, df = 1)#1-pnorm(q=zscores, mean=0,sd=1)#1 - pchisq(q = zscores**2, df = 1)
  
  stat_mu = sum(beta_GH_greta[33:48,z] - beta_GH_greta[33:48,'N2'])
  stat_sd = sqrt(sum(beta_GH_greta_var[33:48,z] + beta_GH_greta_var[33:48,'N2']))
  zscores = stat_mu / stat_sd
  comparisons.hl[3,z] <- 1 - pchisq(q = zscores**2, df = 1)
  
  stat_mu = sum(beta_GH_greta[49:60,z] - beta_GH_greta[49:60,'N2'])
  stat_sd = sqrt(sum(beta_GH_greta_var[49:60,z] + beta_GH_greta_var[49:60,'N2']))
  zscores = stat_mu / stat_sd
  comparisons.hl[4,z] <- 1-pnorm(q=zscores, mean=0,sd=1)#1 - pchisq(q = zscores**2, df = 1)
  
  stat_mu = sum(beta_GH_greta[61:80,z] - beta_GH_greta[61:80,'N2'])
  stat_sd = sqrt(sum(beta_GH_greta_var[61:80,z] + beta_GH_greta_var[61:80,'N2']))
  zscores = stat_mu / stat_sd
  comparisons.hl[5,z] <- 1 - pchisq(q = zscores**2, df = 1)
  
  stat_mu = sum(beta_GH_greta[81:96,z] - beta_GH_greta[81:96,'N2'])
  stat_sd = sqrt(sum(beta_GH_greta_var[81:96,z] + beta_GH_greta_var[81:96,'N2']))
  zscores = stat_mu / stat_sd
  comparisons.hl[6,z] <- 1 - pchisq(q = zscores**2, df = 1)
  
  stat_mu = sum(beta_GH_greta[97:98,z] - beta_GH_greta[97:98,'N2'])
  stat_sd = sqrt(sum(beta_GH_greta_var[97:98,z] + beta_GH_greta_var[97:98,'N2']))
  zscores = stat_mu / stat_sd
  comparisons.hl[7,z] <- 1 - pchisq(q = zscores**2, df = 1)
  
  stat_mu = sum(beta_GH_greta[99:104,z] - beta_GH_greta[99:104,'N2'])
  stat_sd = sqrt(sum(beta_GH_greta_var[99:104,z] + beta_GH_greta_var[99:104,'N2']))
  zscores = stat_mu / stat_sd
  comparisons.hl[8,z] <- 1 - pchisq(q = zscores**2, df = 1)
  
  stat_mu = sum(beta_GH_greta[105:106,z] - beta_GH_greta[105:106,'N2'])
  stat_sd = sqrt(sum(beta_GH_greta_var[105:106,z] + beta_GH_greta_var[105:106,'N2']))
  zscores = stat_mu / stat_sd
  comparisons.hl[9,z] <- 1 - pchisq(q = zscores**2, df = 1)
  
  stat_mu = sum(beta_GH_greta[107:112,z] - beta_GH_greta[107:112,'N2'])
  stat_sd = sqrt(sum(beta_GH_greta_var[107:112,z] + beta_GH_greta_var[107:112,'N2']))
  zscores = stat_mu / stat_sd
  comparisons.hl[10,z] <- 1 - pchisq(q = zscores**2, df = 1)
  
  stat_mu = sum(beta_GH_greta[113:119,z] - beta_GH_greta[113:119,'N2'])
  stat_sd = sqrt(sum(beta_GH_greta_var[113:119,z] + beta_GH_greta_var[113:119,'N2']))
  zscores = stat_mu / stat_sd
  comparisons.hl[11,z] <- 1 - pchisq(q = zscores**2, df = 1)
}

comparisons.hl.adj <- comparisons.hl
for (j in 1:nrow(comparisons.hl)) {
  comparisons.hl.adj[j,] <- p.adjust(comparisons.hl[j,],method='BH')
}
comparisons.hl.adj.all <- matrix(p.adjust(comparisons.hl,method='BH'), nrow = 11, ncol = ncol(beta_GH_greta)-1, 
                                 dimnames = list(rownames(comparisons.hl), colnames(beta_GH_greta)[-1]))


q <- list()
j = 1
thr = 0.05
rate.qvalues <- p.adjust(pv$sub, method = 'BH')
first_points <- c(rep('A*A',6), '2 bp', '1 bp (in repeats)', '1-49 bp', '1 bp (in repeats)', "Tandem dup.")
last_points <- c(rep('T*T',6), 'over 2 bp', "over 50 bp",'50-400 bp',"over 50 bp","Foldbacks")
colors = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE",
           "brown","#8DD3C7","#FFFFB3","#BEBADA","darkmagenta")
names(first_points) = names(last_points) = names(colors) <- c(c("C>A","C>G","C>T","T>A","T>C","T>G"), 'MNV', 'Deletions',
                                                              'Complex\n indels','Insertions', 'Structural\n Variants')
for (proper_name in names(print_names)) {
  
  n <- colnames(count.data)[count.data[data$Sample[data$Genotype.new == proper_name & data$Sample %in% rownames(ma.Y)][1],]>0]
  if (proper_name == "rad(TG3312)-54B")
    n <- "rad.TG3312..54B"
  
  ym = max(beta_GH_greta_high[,n,drop=F])
  
  if (ym < 0.5) ym = 0.5
  
  
  generations <- data$Generation[match(rownames(count.data),data$Sample)][grep(n,CD2Mutant[rownames(count.data)])]
  max_generation <- max(generations)
  print_name <- print_names[proper_name]
  
  if (proper_name == "rad(TG3312)-54B") {
    
    number_of_muts <- 
      rowSums(Y[names(CD2Mutant)[CD2Mutant=='rad-54B:20'],1:119])
    
    number_of_svs <- Y[names(CD2Mutant)[CD2Mutant=='rad-54B:20'],120]
    
    title1 <- paste0('Effects for ',print_name, '; ', round(sum(beta_GH_greta[,n,drop=F]),1), ' (',
                     round(sum(beta_GH_greta_low[,n,drop=F]),1),'-', round(sum(beta_GH_greta_high[,n,drop=F]),1),
                     ') het. mut-s/gen., max 20 gen. with ',
                     round(mean(rowSums(Y)[names(CD2Mutant)[CD2Mutant=='rad-54B:20']])), ' (',
                     min(rowSums(Y)[names(CD2Mutant)[CD2Mutant=='rad-54B:20']]), '-',
                     max(rowSums(Y)[names(CD2Mutant)[CD2Mutant=='rad-54B:20']]), ') mut-s, ',
                     sum(CD2Mutant==paste0('rad-54B',':',max_generation)), #' samples')
                     ' samples, ', round(mean(rowSums(Y[,113:119])[names(CD2Mutant)[CD2Mutant=='rad-54B:20']])), ' (',
                     min(rowSums(Y[,113:119])[names(CD2Mutant)[CD2Mutant=='rad-54B:20']]), '-',
                     max(rowSums(Y[,113:119])[names(CD2Mutant)[CD2Mutant=='rad-54B:20']]), ') SVs per sample')
    
  } else {
    
    number_of_svs <- rowSums(Y[names(CD2Mutant)[CD2Mutant==paste0(proper_name,':',max_generation)],113:119])
    
    number_of_muts <- 
      rowSums(Y[names(CD2Mutant)[CD2Mutant==paste0(proper_name,':',max_generation)],1:119])
    
    
    title1 <- paste0('Effects for ',print_name, '; ', round(sum(beta_GH_greta[,n]),1), ' (',
                     round(sum(beta_GH_greta_low[,n]),1),'-', round(sum(beta_GH_greta_high[,n]),1),
                     ') het. mut-s/gen., max ', max_generation, ' gen. with ',
                     round(mean(number_of_muts)), ' (',
                     min(number_of_muts), '-',
                     max(number_of_muts), ') mut-s, ',
                     sum(CD2Mutant==paste0(proper_name,':',max_generation)), # ' samples')
                     ' samples, ',
                     round(mean(number_of_svs)), ' (',
                     min(number_of_svs), '-',
                     max(number_of_svs), ') SVs per sample')
    
    
  }
  
  if (n == 'N2') {
    significance <- data.frame(N2 = rep(0,119))
  } else {
    significance <- comparisons[,n,drop=F]
  }
  
  plot_a <- plot_fullsig_wb(beta_GH_greta[,n,drop=F], CI = T,
                            low = beta_GH_greta_low[,n,drop=F], ymax = ym,
                            high = beta_GH_greta_high[,n,drop=F], norm = F,
                            colors = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE",
                                       "brown","#8DD3C7","#FFFFB3","#BEBADA","darkmagenta"),
                            significance = significance,
                            rownames = F, size.all=4,
                            err.bar.size = 0.1, thr = thr) + 
    theme(panel.grid = element_blank(), plot.title = element_text(size=10,margin=margin(0,0,0,0)),
          strip.text.y = element_blank(), axis.text.y = element_text(size=10), axis.title.y = element_text(size=10),
          strip.text.x = element_text(size=10),
          axis.ticks = element_line(size=0.5),
          axis.ticks.length = unit(1,"mm")) + 
    ggtitle(title1)
  
  if (ym <= 1) {
    ysegm = 0.1
  } 
  if (ym > 1 & ym <= 2) {
    ysegm = 0.3
  } 
  if (ym > 2 & ym <= 5) {
    ysegm = 0.5
  } 
  if (ym > 5) {
    ysegm = 5
  }
  
  if (n!='N2') {
    
    if (sum(comparisons.hl.adj.all[,n]<thr)>0) {
      df_lines <- data.frame(
        substitution = names(which(comparisons.hl.adj.all[,n]<thr)),
        variable = n,
        x = first_points[names(which(comparisons.hl.adj.all[,n]<thr))], 
        y = -ysegm, 
        xend =last_points[names(which(comparisons.hl.adj.all[,n]<thr))], 
        yend = -ysegm,
        width = 1,
        color = colors[names(which(comparisons.hl.adj.all[,n]<thr))])
      plot_a <- plot_a + geom_segment(data = df_lines, aes(x=x,y=y,yend=yend,xend=xend), col = df_lines$color, size = 0.2)
    }
    
    if (rate.qvalues[n] < thr) {
      df.stars <- data.frame(substitution = 'C>A', variable = n, x = "C*T",
                             y = ym/2, width = 1, value = ym/2)
      plot_a <- plot_a + geom_text(data = df.stars, aes(x = x, y = y), label = "***", size=3)
    }
  }

  q[[j]] <- plot_a
  
  j = j + 1
  
}

GG_save_pdf(q, 'Mutation_accumulation_signatures_with_significances.pdf', 12,5)

######################################BER############################################################

# BER plots

ber <- c('N2','agt.2',"apn.1","exo.3",'ndx.4','parp.1', 'parp.2', 'tdp.1', 'ung.1')
pretty.names <- c('wild-type', 'agt-2', 'apn-1','exo-3', 'ndx-4', 'parp-1', 'parp-2', 'tdp-1', 'ung-1')
which(rate.qvalues[ber]<0.05)
which(comparisons.hl.adj.all[,ber[-1]]<0.05, arr.ind = T)
df_lines <- data.frame(
  substitution = rownames(which(comparisons.hl.adj.all[,ber[-1]]<0.05, arr.ind = T)),
  variable = pretty.names[-1][which(comparisons.hl.adj.all[,ber[-1]]<0.05, arr.ind = T)[,2]],
  x = first_points[rownames(which(comparisons.hl.adj.all[,ber[-1]]<0.05, arr.ind = T))], 
  y = -0.01, 
  xend =last_points[rownames(which(comparisons.hl.adj.all[,ber[-1]]<0.05, arr.ind = T))], 
  yend = -0.01,
  width = 1,
  color = colors[rownames(which(comparisons.hl.adj.all[,ber[-1]]<0.05, arr.ind = T))])

to.show <- beta_GH_greta[,ber,drop=F]; colnames(to.show) <- pretty.names
to.show.low <- beta_GH_greta_low[,ber,drop=F]; colnames(to.show.low) <- pretty.names
to.show.high <- beta_GH_greta_high[,ber,drop=F]; colnames(to.show.high) <- pretty.names
to.show.sign <- cbind(rep(0,119),comparisons[,ber[-1],drop=F]); colnames(to.show.sign) <- pretty.names

df_stars <- data.frame(substitution = 'C>A', variable = unique(pretty.names[c(which(rate.qvalues[ber]<0.05),
                                                                              which(p.adjust(pv$ind, method = 'BH')[ber]<0.05),
                                                                              which(p.adjust(pv$sv, method = 'BH')[ber]<0.05))]), x = "C*T",
                       y = max(to.show.high)/2, width = 1, value = max(to.show.high)/2)

pdf('~/Mutation accumulation/signatures/BER_signatures.pdf', 10, 12)
plot_fullsig_wb(to.show, CI = T,
                low = to.show.low, ymax = max(to.show.high),
                high = to.show.high, norm = F,
                colors = colors, ytitle = 'Number of mutations per generation',
                significance = to.show.sign,
                thr=0.05, err.bar.size = 0.1) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size=10), strip.text.x = element_text(size=12), 
        axis.text.x = element_text(size = 8), axis.title.y = element_text(size = 12)) + 
  geom_segment(data = df_lines, aes(x=x,y=y,yend=yend,xend=xend), col = df_lines$color, size = 1) #+
#geom_text(data = df_stars, aes(x = x, y = y), label = "***", size=6)
dev.off()

######################################NER############################################################

# NER plots

ner <- c('N2','xpa.1',"xpf.1","xpg.1",'xpc.1','csb.1')
pretty.names <- c('wild-type','xpa-1',"xpf-1","xpg-1",'xpc-1','csb-1')
which(rate.qvalues[ner]<0.05)
which(comparisons.hl.adj.all[,ner[-1]]<0.05, arr.ind = T)
df_lines <- data.frame(
  substitution = rownames(which(comparisons.hl.adj.all[,ner[-1]]<0.05, arr.ind = T)),
  variable = pretty.names[-1][which(comparisons.hl.adj.all[,ner[-1]]<0.05, arr.ind = T)[,2]],
  x = first_points[rownames(which(comparisons.hl.adj.all[,ner[-1]]<0.05, arr.ind = T))], 
  y = -0.005, 
  xend =last_points[rownames(which(comparisons.hl.adj.all[,ner[-1]]<0.05, arr.ind = T))], 
  yend = -0.005,
  width = 1,
  color = colors[rownames(which(comparisons.hl.adj.all[,ner[-1]]<0.05, arr.ind = T))])

to.show <- beta_GH_greta[,ner,drop=F]; colnames(to.show) <- pretty.names
to.show.low <- beta_GH_greta_low[,ner,drop=F]; colnames(to.show.low) <- pretty.names
to.show.high <- beta_GH_greta_high[,ner,drop=F]; colnames(to.show.high) <- pretty.names
to.show.sign <- cbind(rep(0,119),comparisons[,ner[-1],drop=F]); colnames(to.show.sign) <- pretty.names

df_stars <- data.frame(substitution = 'C>A', variable = pretty.names[which(rate.qvalues[ner]<0.05)], x = "C*T",
                       y = max(to.show.high)/2, width = 1, value = max(to.show.high)/2)

pdf('~/Mutation accumulation/signatures/NER_signatures.pdf', 10, 10)
plot_fullsig_wb(to.show, CI = T,
                low = to.show.low, ymax = max(to.show.high),
                high = to.show.high, norm = F,
                colors = colors, ytitle = 'Number of mutations per generation',
                significance = to.show.sign,
                thr=0.05, err.bar.size = 0.1) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size=10), strip.text.x = element_text(size=12), 
        axis.text.x = element_text(size = 8), axis.title.y = element_text(size = 12)) + 
  geom_segment(data = df_lines, aes(x=x,y=y,yend=yend,xend=xend), col = df_lines$color, size = 1) +
  geom_text(data = df_stars, aes(x = x, y = y), label = "***", size=6)
dev.off()

#################################Figure3#####################################


repair <- c('N2','agt.2','ung.1','xpa.1', 'xpf.1', 'xpg.1', 'xpc.1')
pretty.names <- c('N2', 'agt-2','ung-1','xpa-1', 'xpf-1', 'xpg-1', 'xpc-1')
df_lines <- data.frame(
  substitution = rownames(which(comparisons.hl.adj.all[,repair[-1]]<0.05, arr.ind = T)),
  variable = pretty.names[-1][which(comparisons.hl.adj.all[,repair[-1]]<0.05, arr.ind = T)[,2]],
  x = first_points[rownames(which(comparisons.hl.adj.all[,repair[-1]]<0.05, arr.ind = T))], 
  y = -0.005, 
  xend =last_points[rownames(which(comparisons.hl.adj.all[,repair[-1]]<0.05, arr.ind = T))], 
  yend = -0.005,
  width = 1,
  color = colors[rownames(which(comparisons.hl.adj.all[,repair[-1]]<0.05, arr.ind = T))])

to.show <- beta_GH_greta[,repair,drop=F]; colnames(to.show) <- pretty.names
to.show.low <- beta_GH_greta_low[,repair,drop=F]; colnames(to.show.low) <- pretty.names
to.show.high <- beta_GH_greta_high[,repair,drop=F]; colnames(to.show.high) <- pretty.names
to.show.sign <- cbind(rep(0,119),comparisons[,repair[-1],drop=F]); colnames(to.show.sign) <- pretty.names

df_stars <- data.frame(substitution = 'C>A', variable = pretty.names[which(rate.qvalues[repair]<0.05)], x = "C*T",
                       y = max(to.show.high)/2, width = 1, value = max(to.show.high)/2)

pdf('BER_NER_selected_signatures.pdf', 10, 10)
plot_fullsig_wb(to.show, CI = T,
                low = to.show.low, ymax = max(to.show.high),
                high = to.show.high, norm = F,
                colors = colors, ytitle = 'Number of mutations per generation',
                significance = to.show.sign,
                thr=0.05, err.bar.size = 0.05) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size=10), strip.text.x = element_text(size=12), 
        axis.text.x = element_text(size = 8), axis.title.y = element_text(size = 12)) + 
  geom_segment(data = df_lines, aes(x=x,y=y,yend=yend,xend=xend), col = df_lines$color, size = 1) +
  geom_text(data = df_stars, aes(x = x, y = y), label = "***", size=6)
dev.off()

#####################################TLS#######################################

tls <- c('N2','polh.1','polh.lf31..1',"rev.gk147834..1","rev.1",'rev.3','polk.1','polq.1')
pretty.names <- c('wild-type','polh-1\n tm3317','polh-1\n (lf31)', 'rev-1\n (gk147834)', 'rev-1\n gk924750', 'rev-3', 'polk-1', 'polq-1')
df_lines <- data.frame(
  substitution = rownames(which(comparisons.hl.adj.all[,tls[-1]]<0.05, arr.ind = T)),
  variable = pretty.names[-1][which(comparisons.hl.adj.all[,tls[-1]]<0.05, arr.ind = T)[,2]],
  x = first_points[rownames(which(comparisons.hl.adj.all[,tls[-1]]<0.05, arr.ind = T))], 
  y = -0.01, 
  xend =last_points[rownames(which(comparisons.hl.adj.all[,tls[-1]]<0.05, arr.ind = T))], 
  yend = -0.01,
  width = 1,
  color = colors[rownames(which(comparisons.hl.adj.all[,tls[-1]]<0.05, arr.ind = T))])

to.show <- beta_GH_greta[,tls,drop=F]; colnames(to.show) <- pretty.names
to.show.low <- beta_GH_greta_low[,tls,drop=F]; colnames(to.show.low) <- pretty.names
to.show.high <- beta_GH_greta_high[,tls,drop=F]; colnames(to.show.high) <- pretty.names
to.show.sign <- cbind(rep(0,119),comparisons[,tls[-1],drop=F]); colnames(to.show.sign) <- pretty.names

df_stars <- data.frame(substitution = 'C>A', variable = pretty.names[c(which(rate.qvalues[tls]<0.05),
                                                                       which(p.adjust(pv$ind, method = 'BH')[tls]<0.05),
                                                                       which(p.adjust(pv$sv, method = 'BH')[tls]<0.05))], x = "C*T",
                       y = max(to.show.high)/2, width = 1, value = max(to.show.high)/2)
pdf('~/Mutation accumulation/signatures/tls_signatures.pdf', 10, 10)
plot_fullsig_wb(to.show, CI = T,
                low = to.show.low, ymax = 0.3,
                high = to.show.high, norm = F,
                colors = colors, ytitle = 'Number of mutations per generation',
                significance = to.show.sign,
                thr=0.05, err.bar.size = 0.1, yspace = 'free_y') +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size=10), strip.text.x = element_text(size=12), 
        axis.text.x = element_text(size = 8), axis.title.y = element_text(size = 12)) + 
  geom_segment(data = df_lines, aes(x=x,y=y,yend=yend,xend=xend), col = df_lines$color, size = 1) +
  geom_text(data = df_stars, aes(x = x, y = y), label = "***", size=6)
dev.off()

############################################ICLR################################

iclr <- c('N2','dog.1','fan.1',"fcd.2","fnci.1",'fncm.1','helq.1')
pretty.names <- c('wild-type','dog-1','fan-1',"fcd-2","fnci-1",'fncm-1','helq-1')
df_lines <- data.frame(
  substitution = rownames(which(comparisons.hl.adj.all[,iclr[-1]]<0.05, arr.ind = T)),
  variable = pretty.names[-1][which(comparisons.hl.adj.all[,iclr[-1]]<0.05, arr.ind = T)[,2]],
  x = first_points[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<0.05, arr.ind = T))], 
  y = -0.01, 
  xend =last_points[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<0.05, arr.ind = T))], 
  yend = -0.01,
  width = 1,
  color = colors[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<0.05, arr.ind = T))])

to.show <- beta_GH_greta[,iclr,drop=F]; colnames(to.show) <- pretty.names
to.show.low <- beta_GH_greta_low[,iclr,drop=F]; colnames(to.show.low) <- pretty.names
to.show.high <- beta_GH_greta_high[,iclr,drop=F]; colnames(to.show.high) <- pretty.names
to.show.sign <- cbind(rep(0,119),comparisons[,iclr[-1],drop=F]); colnames(to.show.sign) <- pretty.names

df_stars <- data.frame(substitution = 'C>A', variable = unique(pretty.names[c(which(rate.qvalues[iclr]<0.05),
                                                                              which(p.adjust(pv$ind, method = 'BH')[iclr]<0.05),
                                                                              which(p.adjust(pv$sv, method = 'BH')[iclr]<0.05))]),
                       x = "C*T",
                       y = c(0.5, 0.2, 0.2, 0.2, 0.2), width = 1, value = c(0.5, 0.2, 0.2, 0.2, 0.2))

pdf('~/Mutation accumulation/signatures/ICLR_signatures.pdf', 10, 10)
plot_fullsig_wb(to.show, CI = T,
                low = to.show.low, ymax = max(to.show.high),
                high = to.show.high, norm = F,
                colors = colors, ytitle = 'Number of mutations per generation',
                significance = to.show.sign,
                thr=0.05, err.bar.size = 0.1, diff_scale = T, diff_limits = c(0.5,1,rep(0.5, 5)),
                yspace = 'free_y') +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size=10), strip.text.x = element_text(size=12), 
        axis.text.x = element_text(size = 8), axis.title.y = element_text(size = 12)) + 
  geom_segment(data = df_lines, aes(x=x,y=y,yend=yend,xend=xend), col = df_lines$color, size = 1) #+
#geom_text(data = df_stars, aes(x = x, y = y), label = "***", size=6)
dev.off()

############################################DSBR################################

thr <- 0.05
dsbr <- c('N2','cku.80','lig.4',"polq.1","brc.1",'brd.1','rfs.1','rip.1','rif.1','smc.5','smc.6','mus.81','slx.1',
          'gen.1','lem.3','smg.1','him.6','wrn.1','rcq.5','rtel.1')
pretty.names <- c('N2','cku-80','lig-4',"polq-1","brc-1",'brd-1','rfs-1','rip-1','rif-1','smc-5','smc-6','mus-81','slx-1',
                  'gen-1','lem-3','smg-1','him-6','wrn-1','rcq-5','rtel-1')
df_lines <- data.frame(
  substitution = rownames(which(comparisons.hl.adj.all[,dsbr[-1]]<thr, arr.ind = T)),
  variable = pretty.names[-1][which(comparisons.hl.adj.all[,dsbr[-1]]<thr, arr.ind = T)[,2]],
  x = first_points[rownames(which(comparisons.hl.adj.all[,dsbr[-1]]<thr, arr.ind = T))], 
  y = -0.01, 
  xend =last_points[rownames(which(comparisons.hl.adj.all[,dsbr[-1]]<thr, arr.ind = T))], 
  yend = -0.01,
  width = 1,
  color = colors[rownames(which(comparisons.hl.adj.all[,dsbr[-1]]<thr, arr.ind = T))])

to.show <- beta_GH_greta[,dsbr,drop=F]; colnames(to.show) <- pretty.names
to.show.low <- beta_GH_greta_low[,dsbr,drop=F]; colnames(to.show.low) <- pretty.names
to.show.high <- beta_GH_greta_high[,dsbr,drop=F]; colnames(to.show.high) <- pretty.names
to.show.sign <- cbind(rep(0,119),comparisons[,dsbr[-1],drop=F]); colnames(to.show.sign) <- pretty.names

df_stars <- data.frame(substitution = 'C>A', variable = unique(pretty.names[c(which(rate.qvalues[dsbr]<0.05),
                                                                              which(ind.qvalues[dsbr]<0.05),
                                                                              which(sv.qvalues[dsbr]<0.05))]), 
                       x = "C*T",
                       y = max(to.show.high)/2, width = 1, value = max(to.show.high)/2)

pdf('DSBR_signatures.pdf', 12, 16)
plot_fullsig_wb(to.show, CI = T,
                low = to.show.low, ymax = 0.5,
                high = to.show.high, norm = F,
                colors = colors, ytitle = 'Number of mutations per generation',
                significance = to.show.sign,
                thr=thr, err.bar.size = 0.1, diff_scale = T, 
                diff_limits = rep(0.5,20),
                yspace = 'free_y') +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size=10), strip.text.x = element_text(size=12), 
        axis.text.x = element_text(size = 8), axis.title.y = element_text(size = 12)) + 
  geom_segment(data = df_lines, aes(x=x,y=y,yend=yend,xend=xend), col = df_lines$color, size = 1) #+
#geom_text(data = df_stars, aes(x = x, y = y), label = "***", size=6)
dev.off()

############################fig4#############################

thr = 0.05
iclr <- c('N2','polh.lf31..1','rev.3','polk.1','dog.1','helq.1')
pretty.names <- c('N2','polh-1\n (lf31)','rev-3', 'polk-1','dog-1','helq-1')
df_lines <- data.frame(
  substitution = rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T)),
  variable = pretty.names[-1][which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T)[,2]],
  x = first_points[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T))], 
  y = -0.01, 
  xend =last_points[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T))], 
  yend = -0.01,
  width = 1,
  color = colors[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T))])

to.show <- beta_GH_greta[,iclr,drop=F]; colnames(to.show) <- pretty.names
to.show.low <- beta_GH_greta_low[,iclr,drop=F]; colnames(to.show.low) <- pretty.names
to.show.high <- beta_GH_greta_high[,iclr,drop=F]; colnames(to.show.high) <- pretty.names
to.show.sign <- cbind(rep(0,119),comparisons[,iclr[-1],drop=F]); colnames(to.show.sign) <- pretty.names

df_stars <- data.frame(substitution = 'C>A', variable = unique(pretty.names[c(which(rate.qvalues[iclr]<0.05),
                                                                              which(p.adjust(pv$ind, method = 'BH')[iclr]<0.05),
                                                                              which(p.adjust(pv$sv, method = 'BH')[iclr]<0.05))]), x = "C*T",
                       y = c(0.4, 0.8, 0.4, 0.4, 0.4), width = 1)

pdf('MA_figure_4_signatures.pdf', 10, 6)
plot_fullsig_wb(to.show, CI = T,
                low = to.show.low, ymax = 0.5,
                high = to.show.high, norm = F,
                colors = colors, ytitle = 'Number of mutations per generation',
                significance = to.show.sign,
                thr=thr, err.bar.size = 0.1,
                diff_scale = T, diff_limits = c(0.5,0.5,0.5,0.5,1,0.5)) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size=10), strip.text.x = element_text(size=12), 
        axis.text.x = element_text(size = 8), axis.title.y = element_text(size = 12)) + 
  geom_segment(data = df_lines, aes(x=x,y=y,yend=yend,xend=xend), col = df_lines$color, size = 1) +
  geom_text(data = df_stars, aes(x = x, y = y), label = "***", size=6)
dev.off()

############################################fig5################################

thr <- 0.05
dsbr <- c('N2',"brc.1",'rfs.1','rip.1','smc.5','smc.6','mus.81','slx.1',
          'him.6','rtel.1')
pretty.names <- c('N2',"brc-1",'rfs-1','rip-1','smc-5','smc-6','mus-81','slx-1',
                  'him-6','rtel-1')
df_lines <- data.frame(
  substitution = rownames(which(comparisons.hl.adj.all[,dsbr[-1]]<thr, arr.ind = T)),
  variable = pretty.names[-1][which(comparisons.hl.adj.all[,dsbr[-1]]<thr, arr.ind = T)[,2]],
  x = first_points[rownames(which(comparisons.hl.adj.all[,dsbr[-1]]<thr, arr.ind = T))], 
  y = -0.01, 
  xend =last_points[rownames(which(comparisons.hl.adj.all[,dsbr[-1]]<thr, arr.ind = T))], 
  yend = -0.01,
  width = 1,
  color = colors[rownames(which(comparisons.hl.adj.all[,dsbr[-1]]<thr, arr.ind = T))])

to.show <- beta_GH_greta[,dsbr,drop=F]; colnames(to.show) <- pretty.names
to.show.low <- beta_GH_greta_low[,dsbr,drop=F]; colnames(to.show.low) <- pretty.names
to.show.high <- beta_GH_greta_high[,dsbr,drop=F]; colnames(to.show.high) <- pretty.names
to.show.sign <- cbind(rep(0,119),comparisons[,dsbr[-1],drop=F]); colnames(to.show.sign) <- pretty.names

df_stars <- data.frame(substitution = 'C>A', variable = unique(pretty.names[c(which(rate.qvalues[dsbr]<0.05),
                                                                              which(p.adjust(pv$ind, method = 'BH')[dsbr]<0.05),
                                                                              which(p.adjust(pv$sv, method = 'BH')[dsbr]<0.05))]), x = "C*T",
                       y = c(0.3,0.3,0.3,1.5,0.6,0.3,0.75,0.3,0.3), width = 1, value = max(to.show.high)/2)

pdf('MA_figure_5_signatures.pdf', 12, 12)
plot_fullsig_wb(to.show, CI = T,
                low = to.show.low, ymax = max(to.show.high),
                high = to.show.high, norm = F,
                colors = colors, ytitle = 'Number of mutations per generation',
                significance = to.show.sign,
                thr=thr, err.bar.size = 0.1,
                diff_scale = T, diff_limits = rep(0.4,10)) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size=10), strip.text.x = element_text(size=12), 
        axis.text.x = element_text(size = 8), axis.title.y = element_text(size = 12)) + 
  geom_segment(data = df_lines, aes(x=x,y=y,yend=yend,xend=xend), col = df_lines$color, size = 1) +
  geom_text(data = df_stars, aes(x = x, y = y), label = "***", size=6)
dev.off()

###################################APOPTOSIS##################################

thr = 0.05
iclr <- c('N2','atm.1',"bub.gt2000..3","bub.ok3437..3","san.1", "ced.3",'ced.4','cep.1',"brc.1",
          "brc.1.ced.3","brc.1.cep.1","him.6", "him.6.ced.3","him.6.ced.4","him.6.cep.1")
pretty.names <- c('wild-type','atm-1',"bub-3\n (gt2000)","bub-3\n (ok3437)","san-1", "ced-3",'ced-4','cep-1',
                  'brc-1\n brd-1', 'brc-1 brd-1;\n ced-3','brc-1 brd-1;\n cep-1','him-6', "him-6\n ced-3","him-6;\n ced-4","him-6;\n cep-1")
df_lines <- data.frame(
  substitution = rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T)),
  variable = pretty.names[-1][which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T)[,2]],
  x = first_points[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T))], 
  y = -0.01, 
  xend =last_points[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T))], 
  yend = -0.01,
  width = 1,
  color = colors[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T))])

to.show <- beta_GH_greta[,iclr,drop=F]; colnames(to.show) <- pretty.names
to.show.low <- beta_GH_greta_low[,iclr,drop=F]; colnames(to.show.low) <- pretty.names
to.show.high <- beta_GH_greta_high[,iclr,drop=F]; colnames(to.show.high) <- pretty.names
to.show.sign <- cbind(rep(0,119),comparisons[,iclr[-1],drop=F]); colnames(to.show.sign) <- pretty.names

df_stars <- data.frame(substitution = 'C>A', variable = unique(pretty.names[c(which(rate.qvalues[iclr]<0.05),
                                                                              which(p.adjust(pv$ind, method = 'BH')[iclr]<0.05),
                                                                              which(p.adjust(pv$sv, method = 'BH')[iclr]<0.05))]), x = "C*T",
                       y = 0.15, width = 1)

pdf('~/Mutation accumulation/signatures/Apoptosis_signatures.pdf', 12, 13)
plot_fullsig_wb(to.show, CI = T,
                low = to.show.low, ymax = 0.5,
                high = to.show.high, norm = F,
                colors = colors, ytitle = 'Number of mutations per generation',
                significance = to.show.sign,
                thr=thr, err.bar.size = 0.1,
                diff_scale = T, diff_limits = c(rep(0.3,8),rep(0.6,3),rep(0.3,4)), yspace = 'free_y') +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size=10), strip.text.x = element_text(size=12), 
        axis.text.x = element_text(size = 8), axis.title.y = element_text(size = 12)) + 
  geom_segment(data = df_lines, aes(x=x,y=y,yend=yend,xend=xend), col = df_lines$color, size = 1) #+
#geom_text(data = df_stars, aes(x = x, y = y), label = "***", size=6)
dev.off()

#################################fig6####################################

thr = 0.05
iclr <- c('N2','atm.1','brc.1','ced.3','cep.1',
          "brc.1.ced.3","brc.1.cep.1",'him.6',"him.6.ced.3")
pretty.names <- c('N2','atm-1',"brc-1 brd-1", "ced-3", "cep-1",
                  'brc-1 brd-1;\n ced-3','brc-1 brd-1;\n cep-1',"him-6","him-6 ced-3")
df_lines <- data.frame(
  substitution = rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T)),
  variable = pretty.names[-1][which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T)[,2]],
  x = first_points[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T))], 
  y = -0.01, 
  xend =last_points[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T))], 
  yend = -0.01,
  width = 1,
  color = colors[rownames(which(comparisons.hl.adj.all[,iclr[-1]]<thr, arr.ind = T))])

to.show <- beta_GH_greta[,iclr,drop=F]; colnames(to.show) <- pretty.names
to.show.low <- beta_GH_greta_low[,iclr,drop=F]; colnames(to.show.low) <- pretty.names
to.show.high <- beta_GH_greta_high[,iclr,drop=F]; colnames(to.show.high) <- pretty.names
to.show.sign <- cbind(rep(0,119),comparisons[,iclr[-1],drop=F]); colnames(to.show.sign) <- pretty.names

df_stars <- data.frame(substitution = 'C>A', variable = unique(pretty.names[c(which(rate.qvalues[iclr]<0.05),
                                                                              which(p.adjust(pv$ind, method = 'BH')[iclr]<0.05),
                                                                              which(p.adjust(pv$sv, method = 'BH')[iclr]<0.05))]), x = "C*T",
                       y = c(0.25, 0.5, 0.5, 0.25, 0.25,0.25), width = 1)

pdf('MA_figure_6_signatures.pdf', 10, 10)
plot_fullsig_wb(to.show, CI = T,
                low = to.show.low, ymax = 0.5,
                high = to.show.high, norm = F,
                colors = colors, ytitle = 'Number of mutations per generation',
                significance = to.show.sign,
                thr=thr, err.bar.size = 0.1,
                diff_scale = T, diff_limits = c(rep(0.5,9)), yspace = 'free_y') +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size=10), strip.text.x = element_text(size=12), 
        axis.text.x = element_text(size = 8), axis.title.y = element_text(size = 12)) + 
  geom_segment(data = df_lines, aes(x=x,y=y,yend=yend,xend=xend), col = df_lines$color, size = 1) #+
#geom_text(data = df_stars, aes(x = x, y = y), label = "***", size=6)
dev.off()

##################################################################

# Session information
# (was run on an LSF for a speed up, hence different versions)

# R version 3.5.1 (2018-07-02)
# Matrix products: default

# attached base packages:
# [1] stats4    parallel  stats     graphics  grDevices utils     datasets
# [8] methods   base

# other attached packages:
# [1] reshape2_1.4.3              ggplot2_3.1.0
# [3] VariantAnnotation_1.28.13   Rsamtools_1.34.1
# [5] Biostrings_2.50.2           XVector_0.22.0
# [7] SummarizedExperiment_1.12.0 DelayedArray_0.8.0
# [9] BiocParallel_1.16.6         matrixStats_0.54.0
# [11] Biobase_2.42.0              GenomicRanges_1.34.0
# [13] GenomeInfoDb_1.18.2         IRanges_2.16.0
# [15] S4Vectors_0.20.1            BiocGenerics_0.28.0
# [17] MASS_7.3-50                 bayesplot_1.7.0
# [19] greta_0.3.0.9002

# loaded via a namespace (and not attached):
# [1] httr_1.4.0               bit64_0.9-7              jsonlite_1.6
# [4] assertthat_0.2.1         blob_1.1.1               BSgenome_1.50.0
# [7] GenomeInfoDbData_1.2.0   progress_1.2.0           globals_0.12.4
# [10] pillar_1.3.1             RSQLite_2.1.1            lattice_0.20-35
# [13] glue_1.3.1               reticulate_1.11.1        digest_0.6.18
# [16] colorspace_1.4-1         Matrix_1.2-14            plyr_1.8.4
# [19] XML_3.98-1.19            pkgconfig_2.0.2          biomaRt_2.38.0
# [22] listenv_0.7.0            zlibbioc_1.28.0          purrr_0.3.2
# [25] scales_1.0.0             whisker_0.3-2            tibble_2.1.1
# [28] withr_2.1.2              GenomicFeatures_1.34.8   lazyeval_0.2.2
# [31] magrittr_1.5             crayon_1.3.4             memoise_1.1.0
# [34] future_1.12.0            tools_3.5.1              prettyunits_1.0.2
# [37] hms_0.4.2                stringr_1.4.0            munsell_0.5.0
# [40] AnnotationDbi_1.44.0     compiler_3.5.1           rlang_0.3.3
# [43] grid_3.5.1               RCurl_1.95-4.12          ggridges_0.5.1
# [46] bitops_1.0-6             base64enc_0.1-3          gtable_0.3.0
# [49] codetools_0.2-15         DBI_1.0.0                R6_2.4.0
# [52] GenomicAlignments_1.18.1 tfruns_1.4               dplyr_0.8.0.1
# [55] tensorflow_1.10          rtracklayer_1.42.2       bit_1.1-14
# [58] stringi_1.4.3            Rcpp_1.0.1               tidyselect_0.2.5
# [61] coda_0.19-2