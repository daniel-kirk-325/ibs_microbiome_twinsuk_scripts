library(NetCoMi)
library(ggplot2)
library(dplyr)
setwd("C:/Users/dan_k/OneDrive - King's College London/Year 2/microbiome_ibs/microbiome_laatste/for_manuscript/phage_analysis/phage_data_2026/results/network_results_hpc")

se_net <- readRDS(file = "original_data/se_net.rds")

# 1. Get node (taxon) names from the network object
nodeNames <- rownames(se_net$input$adjaMat1)  

# 2. Get genus
genus <- sub("_.*", "", nodeNames)

# 3. Format labels
labels <- sub("_", "\n", nodeNames)        # replaces only the first underscore
labels <- gsub("_", " ", labels)           # replaces all remaining underscores
labels <- sub("unclassified", '', labels)  # Remove the word "unclassified" for readability
labels[grepl('phage', labels)] <- ''       # Remove phage labels
names(labels) <- nodeNames


# 4. One colour per genus
genusLevels <- unique(genus)
genusCols   <- rainbow(length(genusLevels))  
names(genusCols) <- genusLevels
genusCols['phage'] <- 'black'              # Convert phage to darkgrey
colorVec <- genusCols[genus]
names(colorVec) <- nodeNames

# 5. Increase number of hubs shown on plot (note hubs on the plot are defined differntly to the manuscript - 
# this is just to align with NetCoMi code)

quant <- 0.90  

evc1 <- se_net$centralities$eigenv1   
evc2 <- se_net$centralities$eigenv2

hub_names1 <- names(evc1)[evc1 >= quantile(evc1, quant)]
hub_names2 <- names(evc2)[evc2 >= quantile(evc2, quant)]

se_net$hubs$hubs1 <- hub_names1
se_net$hubs$hubs2 <- hub_names2

# 6. Plot settings
path <- "C:/Users/dan_k/OneDrive - King's College London/Year 2/microbiome_ibs/microbiome_laatste/for_manuscript/writing/gut_submission/gut_round1/high_res_figures_gut_round1/Figure4.png"
png(
  filename = path,
  width = 18000,
  height = 12000,
  #units = "px",
  res = 900
)

thr <- 0.35

par(oma = c(0, 0, 3, 0))
plot(se_net,
     edgeFilter = "threshold",
     sameLayout = TRUE,
     edgeFilterPar = thr,
     nodeColor = "colorVec",
     colorVec = colorVec,
     nodeTransp = 40,
     hubTransp = 70,
     hubBorderCol = "gray80",   
     hubBorderWidth = 1,      
     labels = labels,
     labelLength = 10,
     borderWidth = 1,
     nodeSize = "eigenvector",
     nodeSizeSpread = 2,
     rmSingles = 'inboth',
     repulsion = 0.6,
     labelScale = FALSE,
     cexLabels = 1.1,
     cexTitle = 1.9,
     title1 = "IBS",
     title2 = "Control",
     showTitle = TRUE)
mtext("Bacteria-Bacteria, Bacteria-vOTU Network Comparison", outer = TRUE, cex = 1.6, line = -0.2)
mtext(paste("(|r| >", thr, ")"), outer = TRUE, cex = 1.4, line = -1.5, col = 'grey40')
dev.off()



# Obtain hubs in each network ---------------------------------------------

thresh <- quantile(se_net$centralities$eigenv1, 0.98)
ibs_hubs <- names(se_net$centralities$eigenv1)[se_net$centralities$eigenv1 >= thresh]

thresh <- quantile(se_net$centralities$eigenv2, 0.98)
ctl_hubs <- names(se_net$centralities$eigenv2)[se_net$centralities$eigenv2 >= thresh]



# Alignment of network correlations with host predictions -----------------

#### Load and processing viral metadata containing host predictions

## Load dataset for mapping phage name to psuedo phage ID
viral_naming_map <- read.csv("C:/Users/dan_k/OneDrive - King's College London/Year 2/microbiome_ibs/microbiome_laatste/for_manuscript/phage_analysis/phage_data_2026/Database_independent_approach_revised/viral contigs metadata/phage_contig_map_with_taxonomy_and_life_cycle.csv")

## Load viral metadata (host, taxonomy, replication cycle information)
host_info <- read.csv("C:/Users/dan_k/OneDrive - King's College London/Year 2/microbiome_ibs/microbiome_laatste/for_manuscript/phage_analysis/phage_data_2026/Database_independent_approach_revised/viral contigs metadata/all_Host_prediction_to_genome.csv")

# Add pseudo phage ID to host_info 
host_info$sudo_phage_ID <- viral_naming_map$sudo_phage_ID[match(host_info$Virus, viral_naming_map$phage_ID)]



#### Format host_info to contain formatted host taxonomy

## Add genus 
host_info$host_genus <- sub(".*g__([^;]*).*", "\\1", host_info$Host.taxonomy)

## Replace underscore plus capital letter with blank space (e.g., Blautia_A --> Blautia)
host_info$host_genus <- gsub("_[A-Z][A-Z]?", "", host_info$host_genus)

## Keep relevant columns and drop repeated entries
host_info <- host_info[, c('sudo_phage_ID', 'host_genus')] %>% 
  distinct()


#### Obtain non-zero correlations for phages in either network and 
#### the bacteria they correlate with

## Get adjacency matrices
assom_ibs <- as.data.frame(se_net$input$assoMat1)
assom_ctl <- as.data.frame(se_net$input$assoMat2)

## Convert to long format (upper triangle)
assom_ibs <- assom_ibs %>%
  tibble::rownames_to_column("Var1") %>%
  tidyr::pivot_longer(cols = -Var1, names_to = "Var2", values_to = "value") %>%
  dplyr::filter(Var1 < Var2) %>% 
  mutate(Group = 'IBS')

assom_ctl <- assom_ctl %>%
  tibble::rownames_to_column("Var1") %>%
  tidyr::pivot_longer(cols = -Var1, names_to = "Var2", values_to = "value") %>%
  dplyr::filter(Var1 < Var2) %>% 
  mutate(Group = 'Ctl')


## Retain associations involving a phage and exceeding |0.35| in either matrix

long_assom <- bind_rows(assom_ibs, assom_ctl) %>% 
  filter(
    abs(value) > 0.35,
    grepl('phage', Var1) | grepl('phage', Var2) 
  )
long_assom <- as.data.frame(long_assom)


## Make sure phages are always in second column
to_swap <- grepl('phage', long_assom$Var1) & !grepl('phage', long_assom$Var2)
tempoary_df <- long_assom$Var1[to_swap]
long_assom$Var1[to_swap] <- long_assom$Var2[to_swap]
long_assom$Var2[to_swap] <- tempoary_df


## Add host information to long_assom

hosts <- c()
for (phageID in long_assom$Var2) {
  hosts <- c(hosts, paste(host_info$host_genus[host_info$sudo_phage_ID == phageID], collapse = ','))
}

long_assom$host_genus <- hosts

write.csv(long_assom, "../network_correlation_hosts.csv")
