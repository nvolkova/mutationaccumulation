# mutation rates

Y <- read.csv('nb5/WORM_FILTERING/Spectrum_119_most_recent_FULL_120819.csv', row.names = 1)
Y[is.na(Y)] <- 0 
Y[,113:119] <- 0

load('nb5/WORM_FILTERING/SV_filtering_DELLY_most_recent_2020.RData')

for (w in intersect(rownames(Y),names(delly.filtcounts))) {
  Y[w,113:119] <- delly.filtcounts[[w]]
}
Y[,120] <- rowSums(Y[,113:119])

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

Y <- Y[setdiff(rownames(Y), c('CD0715b','CD0716b','CD0028a','CD0028c','CD0028d',
                              'CD0093a','CD0093c','CD0093d','CD0245e','CD0245c','CD0245a','CD0245d')),]

data <- data[match(intersect(rownames(Y),names(CD2Mutant)), data$Sample),]
CD2Mutant <- CD2Mutant[intersect(rownames(Y),names(CD2Mutant))]
Y <- Y[names(CD2Mutant),1:119]

counts <- rowSums(Y[is.na(data$Mutagen[match(rownames(Y), data$Sample)]),1:96]) 
generations <- sapply(data$Generation[match(names(counts), data$Sample)], generation.function)
names(generations) <- names(counts)
counts <- counts[!is.na(generations)]
generations <- generations[!is.na(generations)]
genotypes <- t(sapply(names(counts), function(x) as.numeric(unique(data$Genotype.new) == data$Genotype.new[match(x,data$Sample)])))
colnames(genotypes) <- unique(data$Genotype.new)
genotypes <- genotypes[,colSums(genotypes * generations)>0]
counts <- counts[rowSums(genotypes * generations) >0 ]
genotypes <- genotypes[rowSums(genotypes * generations)>0, ]
generations <- generations[names(counts)]
count.data <- data.frame(y = counts+1, genotypes * generations)
count.data <- count.data[,-match(c("trt.1","mrt.2","pot.2"),colnames(count.data))]
count.data <- count.data[rowSums(count.data)>0,]

model <- glm(y ~ ., family = stats::poisson(link = 'identity'), data = count.data)
coef(summary(model)) -> matrix_of_params

names(print_names) <- c(sapply(print_names, function(x) unlist(strsplit(x,split = ' '))[1])[1:15],
                        "polh(lf31)-1",'polk-1',
                        'rev(gk147834)-1',
                        sapply(print_names, function(x) unlist(strsplit(x,split = ' '))[1])[19:32],
                        "rad(TG3308)-54B", "rad(TG3312)-54B",
                        sapply(print_names, function(x) unlist(strsplit(x,split = ' '))[1])[35:44],
                        "bub(gt2000)-3","bub(ok3437)-3",
                        sapply(print_names, function(x) unlist(strsplit(x,split = ' '))[1])[c(47:70)])
genotype_labels <- print_names[sapply(rownames(matrix_of_params), function(x) 
  grep(x,names(print_names))[1])]

matrix_of_params <- matrix_of_params[!is.na(genotype_labels),,drop=F]

o <- order(matrix_of_params[,1])
genotype_labels <- print_names[sapply(rownames(matrix_of_params)[o], function(x) 
  grep(x,names(print_names))[1])]

#pdf('Mutation_rates_all_new.pdf', 8, 5)
#par(mar = c(8,4,4,1))
plot(matrix_of_params[o,1], pch = 16, xaxt = 'n', frame = F, ylab = 'Mut-s per gen.',
     las = 2, xlab = '', cex = 0.5, log = 'y', main = 'Total mutation rates')
abline(h = matrix_of_params['N2',1], lty = 2)
axis(side = 1, at = 1:nrow(matrix_of_params), 
     labels = genotype_labels, 
     las = 2, 
     hadj = 1, tick = F,
     cex.axis = 0.6)
arrows(x0 = 1:nrow(matrix_of_params), 
       y0 = matrix_of_params[o,1] - 2*matrix_of_params[o,2],
       y1 = matrix_of_params[o,1] + 2*matrix_of_params[o,2],
       col ='darkgrey', length=0, lwd=1)
pv <- sapply(rownames(matrix_of_params)[-1], function(z) {
  stat_mu = matrix_of_params[z,1] - matrix_of_params['N2',1]
  stat_sd = sqrt(matrix_of_params[z,2]**2 + matrix_of_params['N2',2]**2)
  zscore = stat_mu / stat_sd
  return(1 - pchisq(q = zscore**2, df = 1))
})
names(which(p.adjust(pv,method='BH')<0.05)) -> o2
arrows(x0 = (1:nrow(matrix_of_params))[match(o2, rownames(matrix_of_params)[o])], 
       y0 = matrix_of_params[o2,1] - 2*matrix_of_params[o2,2],
       y1 = matrix_of_params[o2,1] + 2*matrix_of_params[o2,2],
       col ='darkred', length=0, lwd=2)
legend('topleft', legend = 'significantly different from N2 (FDR 5%)', 
       fill = 'darkred', border = F, bty = "n", cex = 0.7)
#dev.off()

#################################################################3

# More proper way of calculating everything

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

Y[,120] <- rowSums(Y[,113:119])
#Y[c('CD0244a','CD0244c','CD0244d','CD0135d','CD0134d','CD0260e','CD0097c'),113:120] <- 0 
#Y[c('CD0260c'),120] <- 2
#Y[c('CD0260c'),113:119] <- c(0,1,0,0,0,1,0)
#Y[c('CD0259c'),120] <- 2
#Y[c('CD0259c'),113:119] <- c(0,0,0,2,0,0,0)
#Y[c('CD0259c'),113:119] <- c(0,0,0,2,0,0,0)

count.data <- data.frame(genotypes * generations)
ma.Y <- Y[rownames(count.data),]
# ma.Y[c('CD0244a','CD0244c','CD0244d','CD0135d','CD0134d','CD0260e','CD0245c','CD0245d'),113:119] <- 0 
#ma.Y[c('CD0244a','CD0244c','CD0244d','CD0135d','CD0134d','CD0260e','CD0245c','CD0245d','CD0097c'),113:120] <- 0 
#ma.Y[c('CD0244a','CD0244c','CD0244d','CD0135d','CD0134d','CD0260e','CD0097c'),113:120] <- 0 
#ma.Y[c('CD0245e'),120] <- 1
#ma.Y[c('CD0245e'),113:119] <- c(0,1,0,0,0,0,0)
#ma.Y[c('CD0260c'),120] <- 2
#ma.Y[c('CD0260c'),113:119] <- c(0,1,0,0,0,1,0)
#ma.Y[c('CD0259c'),120] <- 2
#ma.Y[c('CD0259c'),113:119] <- c(0,0,0,2,0,0,0)
#ma.Y[c('CD0259c'),113:119] <- c(0,0,0,2,0,0,0)
count.data <- count.data[,-match(c('brc.1.ced.4','brc.1.lig.4',
                                   'rad.TG3312..54B', "rad.TG3308..54B", "mus.81.cep.1", "fcd.2.xpf.1"),colnames(count.data))]
count.data <- count.data[rowSums(count.data)>0,]
ma.Y <- ma.Y[rownames(count.data),]
save(ma.Y, count.data, file = 'yoda1/count.data.full.new.RData')


# on cluster
load('yoda1/count.data.full.RData')
m <- ncol(count.data)
counts <- rowSums(ma.Y)
#library(greta)
sigma <- variable(lower = 0)
rates <-  lognormal(meanlog = 0, sdlog = sigma, dim = c(1, m))
mu = (count.data %*% t(rates))
size = 100
prob = size / (size + mu)
counts <- t(t(counts))
distribution(counts) = negative_binomial(prob = prob, 
                                         size = size * matrix(1,nrow = nrow(count.data), ncol = 1))
ma.model <- model(rates,sigma)
ma.draws <- mcmc(ma.model, warmup = 1000, n_samples = 1000) # do on cluster

library(bayesplot)
mcmc_trace(ma.draws[,grep('rates',colnames(ma.draws[[1]]))[sample(1:m,16)]])
mcmc_trace(ma.draws[,grep('sigma',colnames(ma.draws[[1]])),drop=F])

draws_all <- do.call('rbind',ma.draws)
beta_GH_greta <- matrix(colMeans(draws_all[,grep('rates',colnames(draws_all))]),nrow=ncol(ma.Y)-1)
beta_GH_greta_low <- matrix(apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.025),nrow=ncol(ma.Y)-1)
beta_GH_greta_high <- matrix(apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.975),nrow=ncol(ma.Y)-1)
beta_GH_greta_var <- matrix(apply(draws_all[,grep('rates',colnames(draws_all))], 2, var),nrow=ncol(ma.Y)-1)
colnames(beta_GH_greta) = colnames(beta_GH_greta_high) = colnames(beta_GH_greta_low) = colnames(beta_GH_greta_var) <- colnames(count.data)
save(beta_GH_greta,beta_GH_greta_low,beta_GH_greta_high,beta_GH_greta_var, file = 'full_beta_GH_for_MA_120819.RData')

rates_est<- colMeans(draws_all[,grep('rates',colnames(draws_all))])
rates_low <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.025)
rates_high <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, quantile, 0.975)
rates_var <- apply(draws_all[,grep('rates',colnames(draws_all))], 2, var)

write.csv(rates_est, 'sub_mean_estimates_rates.csv')
write.csv(rates_low, 'sub_low_estimates_rates.csv')
write.csv(rates_high, 'sub_high_estimates_rates.csv')
write.csv(rates_var, 'sub_var_estimates_rates.csv')

# upload
rates_all <- list()
rates_all$est <- read.csv('~/yoda1/mean_estimates_rates_full.csv', row.names = 1)
rates_all$high <- read.csv('~/yoda1/high_estimates_rates_full.csv', row.names = 1)
rates_all$low <- read.csv('~/yoda1/low_estimates_rates_full.csv', row.names = 1)
rates_all$var <- read.csv('~/yoda1/var_estimates_rates_full.csv', row.names = 1)
rownames(rates_all$est) = rownames(rates_all$low) = rownames(rates_all$high) = rownames(rates_all$var) <- 
  colnames(count.data)

rates_sub <- list()
rates_sub$est <- read.csv('~/yoda1/new_rates_for_paper/sub_mean_estimates_rates_full.csv', row.names = 1)
#rates_sub$est <- rbind(rates_sub$est, 0.58, 0.69, 1.148, 0.75, 1.11, 9)

rates_sub$high <- read.csv('~/yoda1/new_rates_for_paper/sub_high_estimates_rates_full.csv', row.names = 1)
#rates_sub$high <- rbind(rates_sub$high,  0.66, 0.67, 1.496, 0.992, 1.32,9.5)

rates_sub$low <- read.csv('~/yoda1/new_rates_for_paper/sub_low_estimates_rates_full.csv', row.names = 1)
#rates_sub$low <- rbind(rates_sub$low, 0.5, 0.51, 0.8, 0.5398, 0.91,8.5)

rates_sub$var <- read.csv('~/yoda1/new_rates_for_paper/sub_var_estimates_rates_full.csv', row.names = 1)
#rates_sub$var <- rbind(rates_sub$var, 0.002, 0.003,0.038, 0.013, 0.01,0.8)

rownames(rates_sub$est) = rownames(rates_sub$low) = rownames(rates_sub$high) = rownames(rates_sub$var) <- 
  colnames(count.data)#, #c('brc-1,lig-4', 'brc-1,ced-4', 
                          #  'rad(gt3308)-54B', 'rad(gt3312)-54B', 
                          #  'fcd-2,xpf-1', 'mus-81,cep-1'))

rates_ind <- list()
rates_ind$est <- read.csv('~/yoda1/new_rates_for_paper/ind_mean_estimates_rates_full.csv', row.names = 1)
#rates_ind$est <- rbind(rates_ind$est, 0.257, 0.15, 0.18, 0.205, 0.185, 0.8)

rates_ind$high <- read.csv('~/yoda1/new_rates_for_paper/ind_high_estimates_rates_full.csv', row.names = 1)
#rates_ind$high <- rbind(rates_ind$high, 0.314, 0.27, 0.3, 0.322, 0.2617, 0.7) 

rates_ind$low <- read.csv('~/yoda1/new_rates_for_paper/ind_low_estimates_rates_full.csv', row.names = 1)
#rates_ind$low <- rbind(rates_ind$low, 0.2, 0.09, 0.06, 0.118, 0.125, 0.9)

rates_ind$var <- read.csv('~/yoda1/new_rates_for_paper/ind_var_estimates_rates_full.csv', row.names = 1)
#rates_ind$var <- rbind(rates_ind$var, 0.001, 0.0018, 0.003, 0.003, 0.0013, 0.0012)

rownames(rates_ind$est) = rownames(rates_ind$low) = rownames(rates_ind$high) = rownames(rates_ind$var) <- colnames(count.data)

rates_sv <- list()
rates_sv$est <- read.csv('~/yoda1/new_rates_for_paper/sv_mean_estimates_rates_full.csv', row.names = 1)
#rates_sv$est <- rbind(rates_sv$est, 0.0385, 0.0386, 0.0629, 0.009, 0.0282, 0.9)

rates_sv$high <- read.csv('~/yoda1/new_rates_for_paper/sv_high_estimates_rates_full.csv', row.names = 1)
#rates_sv$high <- rbind(rates_sv$high, 0.0815, 0.11, 0.131, 0.03, 0.67, 1.02)

rates_sv$low <- read.csv('~/yoda1/new_rates_for_paper/sv_low_estimates_rates_full.csv', row.names = 1)
#rates_sv$low <- rbind(rates_sv$low, 0.0127, 0.001, 0.016, 0.0007, 0.0063, 0.78)

rates_sv$var <- read.csv('~/yoda1/new_rates_for_paper/sv_var_estimates_rates_full.csv', row.names = 1)
#rates_sv$var <- rbind(rates_sv$var, 0.00034, 0.0008, 0.0009, 6.2e-05, 0.00023, 0.005)

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
o <- match(c('N2',DR,BER,NER,TLS,ICLR,DSBR,DS,MMR), rownames(rates_sub$est))
genotype_labels <- print_names[sapply(rownames(rates_sub$est)[o], function(x) 
  grep(x,names(print_names))[1])]

final_colors <- c('#CC6677','#44AA99', '#332288')

pdf('Mutation_rates_6February2020.pdf', 10, 7)
par(mar = c(8,4,4,1))
plot(NA, NA, 
     xlim = c(0,length(genotype_labels)), 
     ylim = log10(c(0.01,215)),
     xaxt = 'n', 
     frame = F, 
     ylab = 'Het. mut-s per gen.',
     las = 2, 
     xlab = '', 
     main = 'Mutation rates', 
     yaxt = 'n')
for (j in seq(2,length(genotype_labels),2))
  polygon(x = c(j-0.3,j-0.3,j+0.7,j+0.7,j-0.3),
          y = log10(c(0.01,215,215,0.01,0.01)),
          col = 'grey86', border = NA)
abline(v = 1.7)
abline(v = 2.7)
abline(v = 9.7)
abline(v = 14.7)
abline(v = 33.7)
abline(v = 40.7)
abline(v = 46.7)
abline(v = 49.7)
abline(v = 58.7)
points(log10(rates_sub$est$x[o]), 
     cex = 0.5, 
     pch = 16, 
     col = final_colors[1])
axis(side = 2, 
     at = log10(c(0.01, 0.1, 0, 0.5, 1, 2, 5, 10, 50, 100, 200)), 
     labels = c('<0.01',0.1, 0, 0.5, 1, 2, 5, 10, 50, 100, 200), 
     las = 2, 
     cex.axis = 0.7)
#polygon(x = c(0,0,63,63,0),
#        y = c(log10(0.5*rates_sub$est['N2','x']),
#              log10(2*rates_sub$est['N2','x']), log10(2*rates_sub$est['N2','x']),
#              log10(0.5*rates_sub$est['N2','x']), log10(0.5*rates_sub$est['N2','x'])),
#        col = rgb(as.vector(col2rgb('darkorange4'))[1]/255,
#                  as.vector(col2rgb('darkorange4'))[2]/255,
#                  as.vector(col2rgb('darkorange4'))[3]/255,0.2), border = NA)
#polygon(x = c(0,0,63,63,0),
#        y = c(log10(0.5*rates_sub$est['N2','x']),
#              log10(2*rates_ind$est['N2','x']), log10(2*rates_ind$est['N2','x']),
#              log10(0.5*rates_ind$est['N2','x']), log10(0.5*rates_ind$est['N2','x'])),
#        col = rgb(as.vector(col2rgb('darkgreen'))[1]/255,
#                  as.vector(col2rgb('darkgreen'))[2]/255,
#                  as.vector(col2rgb('darkgreen'))[3]/255,0.2), border = NA)
#polygon(x = c(0,0,63,63,0),
#        y = c(log10(0.5*rates_sv$est['N2','x']),
#              log10(2*rates_sv$est['N2','x']), log10(2*rates_sv$est['N2','x']),
#              log10(0.5*rates_sv$est['N2','x']), log10(0.5*rates_sv$est['N2','x'])),
#        col = rgb(as.vector(col2rgb('darkslateblue'))[1]/255,
#                  as.vector(col2rgb('darkslateblue'))[2]/255,
#                  as.vector(col2rgb('darkslateblue'))[3]/255,0.2), border = NA)
abline(h = log10(rates_sub$est['N2','x']), lty = 2, col = final_colors[1])
abline(h = log10(rates_ind$est['N2','x']), lty = 2, col = final_colors[2])
abline(h = log10(rates_sv$est['N2','x']), lty = 2, col = final_colors[3])
axis(side = 1, at = 1:nrow(rates_sub$est) + 0.2, 
     labels = genotype_labels, 
     las = 2, 
     hadj = 1, tick = F,
     cex.axis = 0.6, font = 3)
#arrows(x0 = 1:nrow(rates_all$est), 
#       y0 = rates_all$low$x[o],
#       y1 = rates_all$high$x[o],
#       col ='darkgrey', length=0, lwd=1)
arrows(x0 = 1:nrow(rates_sub$est) , 
       y0 = log10(rates_sub$low$x[o]),
       y1 = log10(rates_sub$high$x[o]),
       col =final_colors[1], length=0, lwd=1)
arrows(x0 = 1:nrow(rates_ind$est) + 0.2, 
       y0 = log10(rates_ind$low$x[o]),
       y1 = log10(rates_ind$high$x[o]),
       col =final_colors[2], length=0, lwd=1)
arrows(x0 = 1:nrow(rates_sv$est) + 0.4, 
       y0 = log10(rates_sv$low$x[o]),
       y1 = log10(rates_sv$high$x[o]),
       col =final_colors[3], length=0, lwd=1)
#points(x = c(1:62), y = rates_all$est$x[o], pch = 16, cex = 0.5)
points(x = c(1:62), y = log10(rates_sub$est$x[o]), pch = 16, cex = 0.5, col= final_colors[1])
points(x = c(1:62) + 0.2, y = log10(rates_ind$est$x[o]), pch = 16, cex = 0.5, col = final_colors[2])
yvalues <- rates_sv$est$x[o]; yvalues[yvalues < 0.01] <- 0.01
points(x = c(1:62) + 0.4, y = log10(yvalues), pch = 16, cex = 0.5, col = final_colors[3])
pv <- list()
#pv$all <- sapply(rownames(rates_all$est)[-1], function(z) {
#  stat_mu = rates_all$est[z,'x'] - rates_all$est['N2','x']
#  stat_sd = sqrt(rates_all$var[z,'x'] + rates_all$var['N2','x'])
#  zscore = stat_mu / stat_sd
#  return(1 - pchisq(q = zscore**2, df = 1))
#})
pv$sub <- sapply(rownames(rates_sub$est)[-1], function(z) {
  stat_mu = rates_sub$est[z,'x'] - rates_sub$est['N2','x']
  stat_sd = sqrt(rates_sub$var[z,'x'] + rates_sub$var['N2','x'])
  zscore = stat_mu / stat_sd
  return(1 - pchisq(q = zscore**2, df = 1))
})
pv$ind <- sapply(rownames(rates_ind$est)[-1], function(z) {
  stat_mu = rates_ind$est[z,'x'] - rates_ind$est['N2','x']
  stat_sd = sqrt(rates_ind$var[z,'x'] + rates_ind$var['N2','x'])
  zscore = stat_mu / stat_sd
  return(1 - pchisq(q = zscore**2, df = 1))
})
pv$sv <- sapply(rownames(rates_sv$est)[-1], function(z) {
  stat_mu = rates_sv$est[z,'x'] - rates_sv$est['N2','x']
  stat_sd = sqrt(rates_sv$var[z,'x'] + rates_sv$var['N2','x'])
  zscore = stat_mu / stat_sd
  return(1 - pchisq(q = zscore**2, df = 1))
})
#names(which(p.adjust(pv$sub,method='BH')<0.05)) -> o2
#points(x = (1:nrow(rates_sub$est))[match(o2, rownames(rates_sub$est)[o])],
#       y = rates_sub$est[o2,'x'], col = 'darkred')
#arrows(x0 = (1:nrow(rates_all$est))[match(o2, rownames(rates_all$est)[o])], 
#       y0 = rates_all$low[o2,'x'],
#       y1 = rates_all$high[o2,'x'],
#       col ='darkred', length=0, lwd=2)
names(which(p.adjust(pv$sub,method='BH')<0.05)) -> o2
points(x = (1:nrow(rates_sub$est))[match(o2, rownames(rates_sub$est)[o])] ,
       y = log10(rates_sub$est[o2,'x']), col = c(adjustcolor('darkred', alpha = 0.4),
                                                 'darkred')[as.numeric(rates_sub$est[o2,'x'] > 
                                                                         2 * rates_sub$est['N2','x'] |
                                                                         rates_sub$est[o2,'x'] < 
                                                                         0.5 * rates_sub$est['N2','x']) + 1], 
       pch = 16, cex = 1)
names(which(p.adjust(pv$ind,method='BH')<0.05)) -> o2
points(x = (1:nrow(rates_ind$est))[match(o2, rownames(rates_ind$est)[o])] + 0.2,
       y = log10(rates_ind$est[o2,'x']), col = c(adjustcolor('darkgreen', alpha = 0.4),
                                                 'darkgreen')[as.numeric(rates_ind$est[o2,'x'] > 
                                                                         2 * rates_ind$est['N2','x'] |
                                                                         rates_ind$est[o2,'x'] < 
                                                                         0.5 * rates_ind$est['N2','x']) + 1], 
       pch = 16, cex = 1)
names(which(p.adjust(pv$sv,method='BH')<0.05)) -> o2
points(x = (1:nrow(rates_sv$est))[match(o2, rownames(rates_sv$est)[o])] + 0.4,
       y = log10(rates_sv$est[o2,'x']), col = c(adjustcolor('darkblue', alpha = 0.4),
                                                'darkblue')[as.numeric(rates_sv$est[o2,'x'] > 
                                                                        2 * rates_sv$est['N2','x'] |
                                                                         rates_sv$est[o2,'x'] < 
                                                                        0.5 * rates_sv$est['N2','x']) + 1], 
       pch = 16, cex = 1)
legend('topleft', legend = c('base substitutions','indels','structural variants',
                             'significantly different from N2 (FDR 5%)',
                             'significantly different from N2 (FDR 5%)',
                             'significantly different from N2 (FDR 5%)',
                             'significantly different from N2 (FDR 5%) and over 2-fold difference',
                             'significantly different from N2 (FDR 5%) and over 2-fold difference',
                             'significantly different from N2 (FDR 5%) and over 2-fold difference'), 
       ncol = 1, 
       col = c('darkorange4','darkgreen','darkslateblue',
               adjustcolor('darkred', alpha = 0.4),adjustcolor('darkgreen', alpha = 0.4),adjustcolor('darkblue', alpha = 0.4),
               'darkred','darkgreen','darkblue'), 
       border = F, bty = "n", pch = 16, cex = 0.7, pt.cex = c(rep(0.5,3),rep(1,6)))
dev.off()

# test if the rates are significantly different from N2

# for total rates
p.adjust(pv$sub, method = 'BH') -> rate.qvalues
p.adjust(pv$ind, method = 'BH') -> ind.qvalues
p.adjust(pv$sv, method = 'BH') -> sv.qvalues

#######################################

# cluster frequency
clustdf <- read.csv('summary_of_clusters_per_sample.csv')
clustdf <- clustdf[!(clustdf$genotype %in% c('trt-1','pot-2','mrt-2')),]

clustdf$total_number_of_mutations <- rowSums(Y)[match(clustdf$name, rownames(Y))]

clust <- sapply(rownames(count.data), function(x) {
  if (x %in% clustdf$name) return(clustdf$number_of_clusters[clustdf$name == x])
  else return(0)
})

cd <- count.data[,!(colnames(count.data) %in% c('trt-1','pot-2','mrt-2'))]
cd <- cd[rowSums(cd)>0,]
m <- ncol(cd)
# very dumb model
coef(summary(glm(clust ~ 0 + ., family = stats::poisson(link='identity'), data = cd,
                 start = as.vector(c3)))) -> c1

c3 <- nmSolve(D = matrix(clust, ncol = 1), P = as.matrix(cd))


# Run a model on the number of clusters vs genotype and generation

sigma <- variable(lower = 0)
rates <-  lognormal(meanlog = 0, sdlog = sigma, dim = c(1, m))
mu = (cd %*% t(rates))
clust <- t(t(clust))
distribution(clust) = poisson(lambda = mu)
cl.model <- model(rates,sigma)
cl.draws <- mcmc(cl.model, warmup = 500, n_samples = 500) # do on cluster
draws_all <- do.call('rbind',cl.draws)
rates <- colMeans(draws_all[,1:74])
rates.sd <- apply(draws_all[,1:74],2,sd)
names(rates) = names(rates.sd) <- colnames(cd)

rates <- rates[-c(8,9,32:34,73)]
rates.sd <- rates.sd[-c(8,9,32:34,73)]

pvclust <- NULL
for (zzz in names(rates)[-1]) {
  stat_mu = rates[zzz] - rates['N2']
  stat_sd = sqrt(rates.sd[zzz]**2 + rates.sd['N2']**2)
  zscore = stat_mu / stat_sd
  pvclust <- c(pvclust, 1 - pchisq(q = zscore**2, df = 1))
}
which(p.adjust(pvclust,method='BH') < 0.05)
# him.6, rip.1, brc.1.cep.1, mus.81.cep.1 


################################################################################################

# try testing proportions?
prop <- NULL
N2.total <- sum(ma.Y[CD2Mutant[rownames(ma.Y)] == 'N2:10',1:119])
N2.clust <- sum(clustdf$number_of_clustered_mutations[clustdf$experiment == 'N2:10'])
prop <- c(prop, N2.clust / N2.total)

pv <- NULL
for (zzz in unique(as.character(clustdf$experiment)[-1])) {
  zzz.total <- sum(ma.Y[CD2Mutant[rownames(ma.Y)] == zzz,1:119])
  zzz.clust <- sum(clustdf$number_of_clustered_mutations[clustdf$experiment == zzz])
  pv <- c(pv, fisher.test(rbind(c(N2.clust, N2.total-N2.clust), c(zzz.clust, zzz.total-zzz.clust)), alternative = 'less')$p.value)
  prop <- c(prop, zzz.clust / zzz.total)
  #pv <- c(pv, chisq.test(rbind(c(N2.clust, N2.total-N2.clust), c(zzz.clust, zzz.total-zzz.clust)))$p.value)
}
# fisher: brc-1:40, him-6:20, brc-1,cep-1:20, him-6,ced-4:20, mus-81,cep-1:15, mus-81,cep-1:20
# chi square: him-6:20, brc-1,cep-1:20, him-6,ced-4:20, mus-81,cep-1:15, mus-81,cep-1:20, mlh-1:20, pms-2:20

o <- order(prop)
f <- barplot(prop[o], col = c('white','darkred')[c(1,(p.adjust(pv,method='BH') < 0.05) + 1)][o],
             #ylab = 'Proportion of clustered mutations', las=2,
             cex.axis = 0.7, horiz = T, xaxt = 'n',
             main = 'Proportion of clustered mutations per genotype')
names <- c('N2:10', as.character(unique(clustdf$experiment)[-1]))
axis(side = 2, at = f, labels = names[o], las = 2, cex.axis = 0.7, lwd = 0, lwd.ticks = 1)
axis(side = 3, at = c(0:6) * 0.1, labels = c(0:6) * 0.1, cex.axis = 0.7)
legend('bottomright', legend = 'significantly different from N2 (FDR 5%)', fill = 'darkred',bty = 'n',
       border = NA)
abline(v = prop[1], lty = 2)

#################################

# Run a model on proportions (without generations)

prop.of.clust <- sapply(rownames(cd), function(x) {
  if (x %in% clustdf$name) return(clustdf$number_of_clustered_mutations[clustdf$name == x] / sum(Y[x,c(1:112)]))
  else return(0)
})

sigma <- variable(lower = 0)
sigma2 <- variable(lower = 0)
prop.rates <-  lognormal(meanlog = 0, sdlog = sigma, dim = c(1, m))
cd1 <- cd / rowSums(cd)
mu = (cd1 %*% t(prop.rates))
prop.of.clust <- t(t(prop.of.clust))
distribution(prop.of.clust) = normal(mean = mu, sd = sigma2)
prop.model <- model(prop.rates,sigma,sigma2)
prop.draws <- mcmc(prop.model, warmup = 500, n_samples = 500) # do on cluster
draws_all <- do.call('rbind',prop.draws)
prop.rates <- colMeans(draws_all[,1:74])
prop.rates.sd <- apply(draws_all[,1:74],2,sd)
names(prop.rates) = names(prop.rates.sd) <- colnames(cd)

prop.rates <- prop.rates[-c(8,9,32:34,73)]
prop.rates.sd <- prop.rates.sd[-c(8,9,32:34,73)]

pvprop <- NULL
for (zzz in names(prop.rates)[-1]) {
  stat_mu = prop.rates[zzz] - prop.rates['N2']
  stat_sd = sqrt(prop.rates.sd[zzz]**2 + prop.rates.sd['N2']**2)
  zscore = stat_mu / stat_sd
  pvprop <- c(pvprop, 1 - pchisq(q = zscore**2, df = 1))
}
which(p.adjust(pvprop,method='BH') < 0.05)
# him.6, rip.1, brc.1, agt.2, rfs.1, brc.1.cep.1, mus.81.cep.1

#################################################################

# Visualize

# for a main figure
pdf('~/Cluster_comparison.pdf',16,10)
par(mar = c(2,8,6,2), mfrow = c(1,2))
o <- order(rates)
f <- barplot(rates[o], col = c('white','darkred')[c(1,as.numeric(p.adjust(pvclust,method='BH') < 0.05)+1)][o],
             xaxt = 'n', las = 2, ylab = '', cex.axis = 0.7, yaxt = 'n',
             main = 'Number of clusters per generation', xlim  = c(0,0.5), horiz = T)
points(x = clust[,1] / data$Generation[match(rownames(clust),data$Sample)] + rnorm(nrow(clust),mean=0.01,sd = 0.0001),
       y = f[match(sapply(rownames(clust), function(x) colnames(cd)[cd[x,]>0]), names(rates)[o]),1],
       pch = 16, col = 'gray10', cex = 0.5)
arrows(y0 = f, x0 = rates[o] - 1.96*rates.sd[o], x1 = rates[o] + 1.96*rates.sd[o], length = 0, col = 'gray21', lwd = 0.5)
axis(side=2, at = f, labels = print_names[sapply(names(rates)[o], function(x) grep(x,names(print_names))[1])],
     cex.axis = 0.7, las = 2, font = 3)
axis(side = 3, at = c(0,0.1,0.2,0.3,0.4,0.5), labels = c(0.0,0.1,0.2,0.3,0.4,0.5), cex.axis = 0.7)
abline(v = rates['N2'], lty = 2)

o <- order(prop.rates)
f <- barplot(prop.rates[o], col = c('white','darkred')[c(1,as.numeric(p.adjust(pvprop,method='BH') < 0.05)+1)][o],
             xaxt = 'n', las = 2, ylab = '', cex.axis = 0.7, yaxt = 'n',
             main = 'Proportion of clustered mutations', xlim  = c(0,1), horiz = T)
points(x = prop.of.clust[,1] + rnorm(nrow(prop.of.clust),mean=0.01,sd = 0.0001),
       y = f[match(sapply(rownames(prop.of.clust), function(x) colnames(cd1)[cd1[x,]>0]), names(prop.rates)[o]),1],
       pch = 16, col = 'gray10', cex = 0.5)
arrows(y0 = f, x0 = prop.rates[o] - 1.96*prop.rates.sd[o], x1 = prop.rates[o] + 1.96*prop.rates.sd[o], length = 0,
       col = 'gray21', lwd = 0.5)
axis(side=2, at = f, labels = print_names[sapply(names(prop.rates)[o], function(x) grep(x,names(print_names))[1])],
     cex.axis = 0.7, las = 2, font = 3)
axis(side = 3, at = c(0:10)*0.1, labels = c(0:10)*0.1, cex.axis = 0.7)
abline(v = prop.rates['N2'], lty = 2)
legend('bottomright', legend = 'significantly different from N2 (FDR 5%)', fill = 'darkred',bty = 'n',
       border = NA)
dev.off()


# full figure (for supplement)

pdf('~/Cluster_comparison_main.pdf',10,6)
set.seed(123)
par(mfrow = c(1,2))
boxplot(rates, frame = F, outline = F, ylim = c(0,max(rates+1.96*rates.sd)),
        ylab = 'No. of clusters per generation', main = 'Clusters across genotypes')
#newx <- jitter(rep(1, length(rates)), amount = 0.1)
points(x = newx, y = rates, col = 'gray', pch = 16) 
o <- which(p.adjust(pvclust,method='BH') < 0.1)
points(x = newx[o+1], y = rates[o+1], col = 'darkred', pch = 16) 
arrows(x0 = newx[o+1],y0 = rates[o+1] - 1.96*rates.sd[o+1],y1=rates[o+1]+1.96*rates.sd[o+1],
       col = 'gray21',lwd=0.5,length=0)
abline(h = rates['N2'], lty = 2)
text(x = c(newx[o+1][1] - 0.2, newx[o+1][-1] + 0.2), y = rates[o+1], font = 3,
     labels = print_names[sapply(names(rates)[o+1], function(x) grep(x,names(print_names))[1])], cex = 0.7)

boxplot(prop.rates, frame = F, outline = F, ylim = c(0,max(prop.rates+1.96*prop.rates.sd)),
        ylab = 'Prop. of clustered mut-s', main = 'Proportion of clustered mutations\n across genotypes')
#newx2 <- jitter(rep(1, length(prop.rates)), amount = 0.1)
points(x = newx2, y = prop.rates, col = 'gray', pch = 16) 
o2 <- which(p.adjust(pvprop,method='BH') < 0.1)
points(x = newx2[o2+1], y = prop.rates[o2+1], col = 'darkred', pch = 16) 
arrows(x0 = newx2[o2+1],y0 = prop.rates[o2+1] - 1.96*prop.rates.sd[o2+1],
       y1=prop.rates[o2+1]+1.96*prop.rates.sd[o2+1],
       col = 'gray21',lwd=0.5,length=0)
abline(h = prop.rates['N2'], lty = 2)
text(x = c(newx2[o2+1][1] - 0.2, newx2[o2+1][c(2:3)] + 0.22,newx2[o2+1][4] - 0.2,newx2[o2+1][5] - 0.25,newx2[o2+1][6:7] - 0.2),
     y = prop.rates[o2+1], font = 3,
     labels = print_names[sapply(names(prop.rates)[o2+1], function(x) grep(x,names(print_names))[1])],
     cex = 0.7)
legend('topright', legend = 'significantly different\n from N2 (FDR 5%)', fill = 'darkred',bty = 'n',
       border = NA, cex = 0.7)
dev.off()
