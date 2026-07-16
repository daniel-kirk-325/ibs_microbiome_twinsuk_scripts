library(dplyr) # 1.1.4
library(patchwork)


# CODE EXPLANATION --------------------------------------------------------

## Compiles the results obtained from the network analysis on CREATE HPC.

## First, the network statistics on the observed data are loaded. 

## Then, network statistics on each permuted dataset are loaded and the 
## differences are calculated.

## p-values are derived by calculating the proportion of times the observed 
## difference is more extreme than the permuted difference 


# Compile results from individual permutations and get p-values -----------

overall_stats <- read.csv("original_data/network_statistics.csv")[-1]
overall_stats$diffs <- overall_stats$IBS- overall_stats$Controls 

dir_path <- "permutations"
files <- list.files(path = dir_path, full.names = TRUE)
combined_df <- do.call(rbind, lapply(files, read.csv, header = TRUE))[-1]

# Drop Vertex and Edge connectivity since these are not relevant 
overall_stats <- overall_stats[!overall_stats$Statistic %in% c('Vertex connectivity', 'Edge connectivity'), ]
combined_df <- combined_df[!combined_df$Statistic %in% c('Vertex connectivity', 'Edge connectivity'), ]


# Calculate difference per network statistic in permuted results
# N S2 = 109 (permuted IBS group)
combined_df$diffs <- combined_df$S2- combined_df$S1


## Calculate p-values for each statistic
p_values <- list()
n_sims <- 999
for (stat in unique(overall_stats$Statistic)) {
  null_dist <- combined_df$diffs[combined_df$Statistic == stat]
  obs       <- overall_stats$diffs[overall_stats$Statistic == stat] 
  
  if (obs > 0) {
    p <- 1 - sum(null_dist >= obs) / n_sims
  } else if (obs < 0) {
    p <- 1 - sum(null_dist <= obs) / n_sims
  } else {  # obs == 0
    pos_count <- sum(null_dist > 0)
    neg_count <- sum(null_dist < 0)
    
    if (pos_count > neg_count) {
      p <- sum(null_dist <= obs) / n_sims
    } else if (neg_count > pos_count) {
      p <- sum(null_dist >= obs) / n_sims
    } else {
      p <- 1.0  # symmetric or flat around zero
    }
  }
  
  p_values[[stat]] <- p
}

p_df <- data.frame(Stat = overall_stats$Statistic, Observed_IBS = overall_stats$IBS, Observed_control = overall_stats$Controls, 
                   Observed_diff = overall_stats$diffs, p_value = unlist(p_values)) %>% arrange(p_value)
p_df$q_value <- p.adjust(p_df$p_value, 'fdr')

write.csv(p_df, "permutation_network_statistics_results.csv")

## Visualise results
for (i in unique(p_df$Stat)) {
  hist(combined_df$diffs[combined_df$Statistic == i], main = i, xlab = 'Difference')
  abline(v=overall_stats$diffs[overall_stats$Statistic == i], col = 'red')
}








