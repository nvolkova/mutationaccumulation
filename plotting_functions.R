# Plot a base substitution signature with or without CIs
# mut_matrix - a table with 96-long vectors of base substitutions (in absolute numbers or relative proportions) in columns
# CI = FALSE - an indicator whether the confidence intervals should be drawn
# low, high - if CI = TRUE, these arguments will be screened for matrices containing the 95% low and upper bounds of the entries in mut_matrix
# norm = TRUE - an indicator denoting whether the mut_matrix (and low/high) entries should be normalised before plotting
#               if norm = FALSE, the absolute numbers will be plotted 
plot_sig_wb <- function (mut_matrix, colors = NA, norm=T, ymax = NA, ymin = NA, flip=F, ind=F, font=11, size=12,
                         CI = F, low = NULL, high = NULL, DNV = FALSE, diff_scale = F,
                         rownames = T, diff_limits = NULL, thr = 0.05, size.all = 10, err.bar.size = 0.5) # plotting 96-signatures / profiles with no lines with/out CIs
{
  
  if (is.na(colors)) {
    library(RColorBrewer)
    colors <- c("#2EBAED", "#000000", "#DE1C14","#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  
  IND_TRIPLETS = c(
    "A*A", "A*C", "A*G", "A*T",
    "C*A", "C*C", "C*G", "C*T",
    "G*A", "G*C", "G*G", "G*T",
    "T*A", "T*C", "T*G", "T*T")
  
  ylab = 'Number of mutations'
  if (norm) {
    mult <- colSums(mut_matrix)
    mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    ylab = 'Relative contribution'
  }
  if (is.na(ymax)) ymax = max(mut_matrix)
  if (is.na(ymin)) ymin = 0
  
  
  if (diff_scale)
    scales <- 'free'
  else scales <- 'free_x'
  
  TRIPLETS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
  TRIPLETS_112 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3),IND_TRIPLETS)
  if (norm) norm_mut_matrix = apply(mut_matrix, 2, function(x) x/sum(x)) else norm_mut_matrix <- mut_matrix
  if (missing(colors)) {
    if (!ind)  colors = c("#2EBAED", "#000000", "#DE1C14",
                          "#D4D2D2", "#ADCC54", "#F0D0CE")
    else colors = c("#2EBAED", "#000000", "#DE1C14", "yellow",
                    "#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  context = TRIPLETS_96
  substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16)
  if (ind) {
    context = TRIPLETS_112
    substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G","INDEL"), each = 16)
  }
  substring(context, 2, 2) = "*"
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  if (flip==T) norm_mut_matrix = -norm_mut_matrix
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  if (ymax>0) breaks = c(0,0.1) else breaks = c(-0.1,0)
  if (ymax>0) lbls = c('0','0.1') else lbls = c('0.1','0')

  
  if (!is.null(diff_limits)) {
    df3 <- rbind(df3, data.frame(substitution = 'Zfake', context = '0', variable = unique(df3$variable), value = diff_limits))
  }
  
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution)) +
    geom_bar(stat = "identity", colour = "black", size = NA) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution, scales = scales) +
    ylab(ylab) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = size.all,vjust = 1), 
          axis.text.y = element_text(size = size.all), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = size.all/2, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = size.all), 
          strip.text.y = element_text(size = size), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  if (!diff_scale) plot <- plot + coord_cartesian(ylim = c(0,ymax), clip='off') 
  if (!rownames) plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  

  
  if (CI) {
    norm_mut_matrix_lower <- as.matrix(low)
    norm_mut_matrix_upper <- as.matrix(high)
    if (norm)
      for (i in 1:ncol(norm_mut_matrix_lower)) {
        norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
        norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_upper)[i]]
      }
    rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
    df_lower = cbind(df,as.data.frame(norm_mut_matrix_lower))
    df_upper = cbind(df,as.data.frame(norm_mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    df3$value_min = df_lower$value
    df3$value_max = df_upper$value
    
    if (!is.null(diff_limits)) {
      df3$value_min = c(df_lower$value, rep(0,ncol(mut_matrix)))
      df3$value_max = c(df_upper$value, rep(0,ncol(mut_matrix)))
    } else {
      df3$value_min = df_lower$value
      df3$value_max = df_upper$value
    }
    
    plot <- plot + geom_linerange(data=df3, aes(ymin = value_min, ymax = value_max), colour = 'darkgrey', size = err.bar.size)
  }
  
  return(plot)
  
  
}

############################################################################################################

# visualisation for the huge profile!

TRIPLETS = c(
  "A*A", "A*C", "A*G", "A*T",
  "C*A", "C*C", "C*G", "C*T",
  "G*A", "G*C", "G*G", "G*T",
  "T*A", "T*C", "T*G", "T*T")
SVS = c("Tandem dup.","Deletions","Inversions","Complex events","Translocations","Interchr. events","Foldbacks")
INDELS = c(
  '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
  "6-50 bp","over 50 bp",
  '1-49 bp','50-400 bp',
  '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
  '6-50 bp','over 50 bp')
library(grid)
library(ggplot2)
library(reshape2)

# Plot a full mutation signature with or without CIs
# mut_matrix - a table with 119-long vectors of mutation counts/proportions per column
# CI = FALSE - an indicator whether the confidence intervals should be drawn
# low, high - if CI = TRUE, these arguments will be screened for matrices containing the 95% low and upper bounds of the entries in `mut_matrix`
# norm = TRUE - an indicator denoting whether the mut_matrix (and low/high) entries should be normalised before plotting
#               if norm = FALSE, the raw numbers from mut_matrix will be plotted 
# diff_scale = FALSE - allow for different scales in visulasing different columns of `mut_matrix`
# rownames = TRUE - indicator whether the context names should be plotted (on the bottom of the plot)
# diff_limits = NULL - can be supplied as a list of y-axis limits for plotting different entry columns of `mut_matrix`
#                      only works when `diff_scale = TRUE`
# significance = NULL - can be supplied as a matrix of p-values of comparison to the first column in `mut_matrix`
# thr = 0.05 - threshold to consider the difference significant (as per p-values in `significance`)
plot_fullsig_wb <- function(mut_matrix, colors=NA, norm=T, ymax=NA, ytitle = 'Number of mutations', # signatures for 119 profile, DNV collapsed, indel resolved
                               CI = F, low = NA, high = NA, diff_scale = F, rownames = T, diff_limits = NULL,
                               significance = NULL, thr = 0.05, size.all = 10, err.bar.size = 0.5, yspace = NULL)
{
  types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[',
                      rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']',
                      rep(c('A','C','G','T'), 24), sep='')
  
  ID <- c('Deletions','Complex\n indels','Insertions')
  
  if (is.na(colors[1])) {
    library(RColorBrewer)
    colors <- c(brewer.pal(6, "Set1"), 'brown', brewer.pal(3,"Set2"), 'darkmagenta')
  }
  
  TRIPLETS = c(
    "A*A", "A*C", "A*G", "A*T",
    "C*A", "C*C", "C*G", "C*T",
    "G*A", "G*C", "G*G", "G*T",
    "T*A", "T*C", "T*G", "T*T")
  
  INDELS = c(
    '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
    "6-50 bp","over 50 bp",
    '1-49 bp','50-400 bp',
    '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
    '6-50 bp','over 50 bp')
  
  SVS = c("Tandem dup.","Deletions","Inversions","Complex events","Translocations","Interchr. events","Foldbacks")
  
  if (norm) {
    mult <- colSums(mut_matrix)
    norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    ytitle <- 'Relative contribution'
  } else norm_mut_matrix <- mut_matrix
  
  
  if (is.na(ymax)) ymax <- max(norm_mut_matrix)
  if(is.null(yspace)) yspace = 'fixed'
  
  TRIPLETS_96 = rep(TRIPLETS, 6)
  context = c(TRIPLETS_96,'2 bp','over 2 bp',INDELS,SVS)
  substitution = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16),
                   'MNV',
                   'MNV',
                   rep(ID,c(6,2,6)),
                   rep('Structural\n Variants',7))
  df = data.frame(substitution = substitution, context = context)
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  df3$context <- as.character(df3$context)
  df3$context <- factor(df3$context, levels = unique(context))
  if (!is.null(significance)) {
    significance <- significance[,colnames(mut_matrix)]
    dfsig <- melt(significance)
    df3$significance <- dfsig$value
    df3$significance[df3$significance >= thr] <- 0.4
    df3$significance[df3$significance < thr] <- 1
  } else
    df3$significance <- 1
  df3$substitution = factor(df3$substitution, levels=c("C>A","C>G","C>T","T>A","T>C","T>G",'MNV',ID,"Structural\n Variants"))
  df3$width = rep(0.7,nrow(df3));
  df3$width[df3$substitution == 'MNV' | df3$substitution == 'Complex\n indels'] <- 0.2
  if (diff_scale)
    scales <- 'free'
  else scales <- 'free_x'
  
  if (!is.null(diff_limits)) {
    df3 <- rbind(df3, data.frame(substitution = 'Zfake', context = '0', variable = unique(df3$variable), value = diff_limits, significance = 0.4, width = 0.2))
  }
  
  p = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = width)) +
    geom_bar(stat = "identity", colour = "black", size = NA, alpha = df3$significance[order(df3$variable)]) +
    scale_fill_manual(values = c(colors,'white')) + facet_grid(variable ~ substitution, scales = scales, space = yspace) +
    ylab(ytitle) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = size.all,vjust = 1), 
          axis.text.y = element_text(size = size.all), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = size.all/2, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = size.all), 
          strip.text.y = element_text(size = size.all),
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  if (!diff_scale) p <- p + coord_cartesian(ylim = c(0,ymax), clip='off') 
  if (!rownames) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  
  if (CI) {
    norm_mut_matrix_lower <- as.matrix(low)
    norm_mut_matrix_upper <- as.matrix(high)
    if (norm)
      for (i in 1:ncol(norm_mut_matrix_lower)) {
        norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
        norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_upper)[i]]
      }
    rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
    df_lower = cbind(df,as.data.frame(norm_mut_matrix_lower))
    df_upper = cbind(df,as.data.frame(norm_mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    if (!is.null(diff_limits)) {
      df3$value_min = c(df_lower$value, rep(0,ncol(mut_matrix)))
      df3$value_max = c(df_upper$value, rep(0,ncol(mut_matrix)))
    } else {
      df3$value_min = df_lower$value
      df3$value_max = df_upper$value
    }
    p <- p + geom_linerange(data=df3, aes(ymin = value_min, ymax = value_max), colour = 'darkgrey', size = err.bar.size)
#    plot <- plot + geom_pointrange(data=df3, aes(ymin=value, ymax=value_max, colour=substitution),
#                                  fatten = 0.001, size=0.5, show.legend = F) +
#      scale_color_manual(values=c(colors,"white")) +
#      geom_pointrange(data=df3, aes(ymin=value_min,ymax=value),
#                      size=0.5, fatten = 0.001,show.legend = F, colour="white")
  }

  return(p)
  
}
