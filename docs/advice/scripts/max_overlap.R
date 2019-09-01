# max_overlap script for coverage calculation

# Author: Shee, Zhi Qiang
# Conditions of use: GPLv3
# max_overlap - this script calculates a set of three statistics to estimate capture coverage of target sequences rescued by HybPiper (Johnson et al, 2016)
# I used it to identify samples with low coverage that may affect downstream analysis.
# Statistic 1: representedness = proportion of species/genes for which sequences were obtained
# Statistic 2: completeness = proportion of target sequence obtained for each species/gene
# Statistic 3: evenness = how evenly the sequence lengths are distributed across species/genes, adapted from a measure of species evenness (Pielou EC. 1966. The measurement of diversity in different types of biological collections. J Theor Biol 13: 131-144.)
# coverage = (representedness x completeness x evenness)^(1/3) <-- cube root because otherwise the product of three fractions would be an incredibly small number
# Theoretically, I prefer this approach because it provides a better idea of which combinations of species/genes would work, instead of pure taxon occupancy (representedness) or capture success (completeness), which do not capture evenness.
# Empirically, this approach has NOT been tested.
# If you find any bugs, please modify accordingly and let me know (totally optional). I will probably not have any time to troubleshoot in the near future (maybe when I'm finally retired). Good luck! ZQ

# Input file is the tab-delim txt output of the get_seq_lengths.py HybPiper script, which should be located in the working directory unless otherwise specified
# Output file is a csv and will be written to same directory unless otherwise specified.
workdir <- '/home/zhi/zq/Schefflera/HybPiper' # specify working directory
infile <- 'seq_lengths.txt' # specify input file name
outfile <- 'max_overlap.csv' # specify output file name

# ----- modifications below this line typically not required -----

setwd(workdir)

# create two matrices to store results

sl <- read.table(infile,
                 row.names = 1,
                 header = TRUE)

sp <- sl
for(i in colnames(sl)){
  sp[,i] <- sl[,i] / sl['MeanLength',i]
  sp[sp[,i] > 1,i] <- 1
}

sl <- sl[-1,]
sp <- sp[-1,]

sp <- data.frame(sp,
                 rep_sp = numeric(nrow(sp)),
                 com_sp = numeric(nrow(sp)),
                 eve_sp = numeric(nrow(sp)),
                 cov_sp = numeric(nrow(sp))
                 )

sp2 <- data.frame(matrix(ncol = ncol(sp), nrow = 4), row.names = c('rep_gn', 'com_gn', 'eve_gn', 'cov_gn'))
colnames(sp2) <- colnames(sp)

sp <- rbind(sp, sp2)

# calculate statistics per gene

for(i in 1:ncol(sl)){
  sp[nrow(sl) + 1,i] <- sum(sp[1:nrow(sl),i] != 0) / nrow(sl)
  sp[nrow(sl) + 2,i] <- mean(sp[1:nrow(sl),i])
  x <- sp[1:nrow(sl),i]
  x <- x [x != 0]
  sp[nrow(sl) + 3,i] <- - sum((x/sum(x)) * log(x/sum(x))) / log(nrow(sl))
}

sp[nrow(sl) + 4,] <- (sp[nrow(sl) + 1,] * sp[nrow(sl) + 2,] * sp[nrow(sl) + 3,])^(1/3)

# calculate statistics per sample

for(i in 1:nrow(sl)){
  sp[i,ncol(sl) + 1] <- sum(sp[i,1:ncol(sl)] != 0) / ncol(sl)
  sp[i,ncol(sl) + 2] <- mean(as.numeric(sp[i,1:ncol(sl)]))
  x <- sp[i,1:ncol(sl)]
  x <- x [x != 0]
  sp[i,ncol(sl) + 3] <- - sum((x/sum(x)) * log(x/sum(x))) / log(ncol(sl))
}

sp[,ncol(sl) + 4] <- (sp[,ncol(sl) + 1] * sp[,ncol(sl) + 2] * sp[,ncol(sl) + 3])^(1/3)

write.csv(sp, outfile)
