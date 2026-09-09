library(stringr)
library(ggplot2)
library(dplyr)
library(ggtext)
library(tidyr)
library(ggforce)
library(tidyr)
library(tibble)
library(circlize)

# Explanation  ------------------------------------------------------------

# Uses the results generated from the phage analyses to generate figures and 
# perform additional analyses, namely:

#     - Obtain replication cycle information of IBS-associated phages
#     - Identify hosts of IBS-associated phages
#     - Obtain information regarding AMGs in IBS-associated phages


# Load data ---------------------------------------------------------------

## Load dataset for mapping phage name to psuedo phage ID
viral_naming_map <- read.csv("data_deposition/gut_paper/gut_paper_datasets/metagenome_data/phage_data/Database_independent_approach_revised/viral contigs metadata/phage_contig_map_with_taxonomy_and_life_cycle.csv")


## Load phage x IBS LMM results
all_results <- read.csv("data_deposition/gut_paper/gut_paper_results/phage_results/phage_lmm.csv")

# Get results which passed FDR
res_df <- all_results %>%   filter(q_value < 0.05, !is.na(q_value))


## Load viral metadata (host, taxonomy, replication cycle information)
host_info <- read.csv("data_deposition/gut_paper/gut_paper_datasets/metagenome_data/phage_data/Database_independent_approach_revised/viral contigs metadata/all_Host_prediction_to_genome.csv")
taxonomy <- read.csv("data_deposition/gut_paper/gut_paper_datasets/metagenome_data/phage_data/Database_independent_approach_revised/viral contigs metadata/phage_contig_map_with_taxonomy_and_life_cycle.csv")

# Add pseudo phage ID to host_info 
host_info$sudo_phage_ID <- viral_naming_map$sudo_phage_ID[match(host_info$Virus, viral_naming_map$phage_ID)]


#### Format host_info to contain formatted host taxonomy

## Add genus 
host_info$host_genus <- sub(".*g__([^;]*).*", "\\1", host_info$Host.taxonomy)

## Replace underscore plus capital letter with blank space (e.g., Blautia_A --> Blautia)
host_info$host_genus <- gsub("_[A-Z][A-Z]?", "", host_info$host_genus)


## Add family 
host_info$host_family <- sub(".*f__([^;]*).*", "\\1", host_info$Host.taxonomy)


## Make subsets including only IBS-associated phages
hosts_ibs <- host_info[host_info$sudo_phage_ID %in% res_df$Virus, ]
taxonomy_ibs <- taxonomy[taxonomy$sudo_phage_ID %in% res_df$Virus, ]


## Load bacterial dataframe
top_bacteria_df <- read.csv("data_deposition/gut_paper/gut_paper_results/bacterial_results/bacteria_analysis.csv")
top_bacteria <- top_bacteria_df[top_bacteria_df$q_value < 0.05, ]


# Replication cycle information -------------------------------------------

table(taxonomy_ibs$life_cycle)

# Temperate  Virulent 
# 42         8 

# Remove phages with equal probability of being temperate or viral from the taxonomy dataset
taxonomy_filt <- taxonomy[taxonomy$life_cycle != 'EQUAL', ]

# Obtain non-IBS associated phages that were tested in species analysis
non_ibs_assoc_species <- all_results$Virus[all_results$q_value > 0.05]

## Get lifestyles of IBS vs non-IBS-associated phages
IBS_assoc_rep_cycle <- table(taxonomy_ibs$life_cycle)
non_IBS_assoc_rep_cycle <- table(taxonomy_filt$life_cycle[taxonomy_filt$sudo_phage_ID %in% non_ibs_assoc_species])

# Build contingency matrix 
count_matrix <- data.frame(
  non_IBS = as.vector(non_IBS_assoc_rep_cycle[ names(IBS_assoc_rep_cycle) ]),
  IBS = as.vector(IBS_assoc_rep_cycle),
  row.names = names(IBS_assoc_rep_cycle)
)


## Perform Chi-squared test
chi_sq <- chisq.test(count_matrix)
chi_sq$observed - chi_sq$expected 


# Standardized residuals
std_residuals <- chi_sq$stdres
print("Standardized residuals:")
print(std_residuals)

prop_ci <- prop.table(as.matrix(count_matrix), margin=2)  # proportions by column
print("Column-wise proportions:")
print(prop_ci)


# Get odds ratios with CI (for 2x2 tables only)
if(nrow(count_matrix) == 2 & ncol(count_matrix) == 2){
  or_result <- epitools::oddsratio(as.matrix(count_matrix))
  print("Odds ratio with CI:")
  print(or_result)
}

# Format label with OR and p-value for plot
plot_label <- paste0('OR(95%CI) = ', round(or_result$measure[2,1], 2), 
                ' (', round(or_result$measure[2,2],2), '-', round(or_result$measure[2,3],2), ')',
                '\nChi-squared p-value = ', signif(chi_sq$p.value, 2), '\n*')

## OR = 0.45

#### Make plot
rep_cycle_prop_long <- as.data.frame(prop_ci) %>%
  rownames_to_column("rep_cyc") %>%
  pivot_longer(-rep_cyc, names_to = "group", values_to = "proportion")

rep_cycle_prop_long$group <- ifelse(rep_cycle_prop_long$group == 'IBS', 'IBS-associated\n(FDR<0.05)', 'non-IBS-associated\n(FDR>0.05)')

plot <- ggplot(rep_cycle_prop_long, aes(x = group, y = proportion, fill = rep_cyc)) +
  geom_col(position = "fill", color = "black") +
  labs(
    title = "\nReplication Cycle in IBS-associated\nvs non-IBS-associated vOTUs\n",
    x = "",
    y = "Proportion\n",
    fill = ""
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    labels = function(x) ifelse(x <= 1, x, ""),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 1.15), clip = "off") +
  scale_fill_viridis_d(option = "viridis", direction = 1) +
  theme_minimal(base_size = 19) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.text.x = element_text(size = 19),
    axis.text.y = element_text(size = 22),
    axis.title.y = element_text(size = 23),
    plot.margin = margin(1.5, 2, 2, 2),
    panel.spacing = unit(0, "lines")
    ) +
  annotate("segment", x = 1, xend = 2, y = 1.04, yend = 1.04, size = 0.8) +
  annotate("segment", x = 1, xend = 1, y = 1.00, yend = 1.04, size = 0.8) +
  annotate("segment", x = 2, xend = 2, y = 1.00, yend = 1.04, size = 0.8) +
  annotate("text", x = 1.5, y = 1.09, label = plot_label, size = 6)



# Circos plot of IBS-associated phages -----------------------------------------

### 1. Label formatting and ordering 
plot_df <- merge(res_df, hosts_ibs[, c('sudo_phage_ID', 'host_genus')], by.x = 'Virus', by.y = 'sudo_phage_ID')

## Keep only one instance of each phage-host combination
plot_df <- plot_df[!duplicated(plot_df[, c('Virus', 'host_genus')]), ]

## Replace underscores with spaces in phage names
plot_df$Virus <- gsub("_", ' ', plot_df$Virus)

## Establish order
ordering_hosts <- sort(table(plot_df$host_genus))


#### Modify phage names which are repeated
repeated_phage <- names(table(plot_df$Virus)[table(plot_df$Virus) == 1])

plot_df <- plot_df %>%
  group_by(Virus) %>%
  mutate(
    dup_count = n(),
    row_num = row_number()  # order within the group (preserves current order)
  ) %>%
  ungroup() %>%
  mutate(
    Virus_label = ifelse(
      dup_count == 1 | row_num == 1,
      Virus,  # uniques and first duplicate keep original
      paste0(Virus, " (", row_num, ")")  # second gets (2), third (3), etc.
    )
  ) %>%
  select(-dup_count, -row_num)

## Format phage names 
plot_df$Virus_label <- paste0("<i>", plot_df$Virus_label, "</i>")

plot_df <- plot_df %>%
  arrange(host_genus, Virus_label) %>%
  mutate(Virus_label = as.character(Virus_label))

ordering_viruses <- unique(plot_df$Virus_label)  
ordering_hosts <- names(sort(table(plot_df$host_genus)))

sector_order <- c(ordering_viruses, ordering_hosts)

### 2. Links

# Assuming 'color_group' is consistent within each Virus (as it should be)
link_color_df <- plot_df %>%
  group_by(Virus_label, host_genus) %>%
  summarise(link_col = first(ifelse(OR > 1, "#FF000088", "#0000FF88")), .groups = "drop")

df_links <- plot_df %>%
  select(host_genus, Virus_label) %>%
  mutate(value = 1) %>%
  group_by(host_genus, Virus_label) %>%
  summarise(value = n(), .groups = "drop") %>%
  left_join(link_color_df, by = c("Virus_label", "host_genus"))

link_colors <- df_links$link_col  # Now matches exactly


### 3. Sector colors

host_colors <- rainbow(length(ordering_hosts))
names(host_colors) <- ordering_hosts

virus_colors <- rep("grey30", length(ordering_viruses))
names(virus_colors) <- ordering_viruses

grid_col <- c(virus_colors, host_colors)


### 4. Gaps

gap_deg <- rep(0.5, length(sector_order))
gap_deg[length(sector_order)] <- 15
gap_deg[length(ordering_viruses)] <- 15


### 5. Initialize 

png(
  "phage_circos.png",
  width = 4000,
  height = 6000,
  res = 900
)

circos.clear()
circos.par(
  start.degree = 90,
  gap.after = gap_deg,
  track.height = 0.08,
  canvas.xlim = c(-1.5, 1.2),   
  canvas.ylim = c(-1.5, 1.2),
  points.overflow.warning = FALSE
)


###  6. Chord diagram

chordDiagram(
  df_links[, c("Virus_label", "host_genus", "value")],  
  order = sector_order,
  grid.col = grid_col,
  col = link_colors,
  #transparency = 0,  
  annotationTrack = "grid",
  annotationTrackHeight = mm_h(3)
)


###  7. OUTER LABEL TRACK 

circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = mm_h(12),
  bg.border = NA,
  panel.fun = function(x, y) {
    
    sector <- get.cell.meta.data("sector.index")
    xcenter <- get.cell.meta.data("xcenter")
    
    # Hosts (right side)
    if (sector %in% ordering_hosts) {
      circos.text(
        xcenter,
        1.65,                  
        sector,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 0.55,
        font = 2,
        clipping = FALSE      
      )
    }
    
    # Viruses 
    if (sector %in% ordering_viruses) {
      label <- sub("</i>$", "", sub("^<i>", "", sector))
      circos.text(
        xcenter,
        1.55,                  
        label,
        facing = "reverse.clockwise",
        niceFacing = TRUE,
        adj = c(1, 0.5),
        cex = 0.25,
        clipping = FALSE
      )
    }
  }
)

###  8. Title + legend

title("IBS-associated vOTUs (FDR<0.05) \nGrouped by Predicted Host Genus", cex.main = 1, line = -2)

legend(
  x = "bottomleft",
  inset = c(0.35, 0),
  legend = c("OR > 1", "OR < 1"),
  col = c("red", "blue"),
  lty = 1,
  lwd = 2,
  bty = "n",
  cex = 1,
  xpd = TRUE
)

dev.off()
circos.clear()




# Host genus barchart ---------------------------------------------------

## Add host genus and family to OR table
OR_hosts <- merge(res_df, hosts_ibs, by.x = 'Virus', by.y = 'sudo_phage_ID', keep.all = TRUE)

## Keep only one instance of each phage-host combination
OR_hosts <- OR_hosts[!duplicated(OR_hosts[, c('Virus', 'host_genus')]), ]


## Add family to unknown genera 
unknown_genera <- c('CAG', 'UBA', 'SFHK')
OR_hosts[grepl(paste0(unknown_genera, collapse = '|'), OR_hosts$host_genus), 'host_genus'] <- paste0('(', OR_hosts$host_family[grepl(paste0(unknown_genera, collapse = '|'), OR_hosts$host_genus)], ') ', OR_hosts[grepl(paste0(unknown_genera, collapse = '|'), OR_hosts$host_genus), 'host_genus'])


## Add label for direction of change 
OR_hosts <- OR_hosts %>% mutate(OR_Type = ifelse(OR > 1, "Positive", "Negative")) 


#### Get order 

## Make table containing positive and negative associations counts
host_OR_counts <- OR_hosts %>% count(host_genus, OR_Type)

## Establish ordering (enforced after adding text labels)
ordering <- host_OR_counts %>% group_by(host_genus) %>% summarise(sum(n)) %>% arrange(`sum(n)`) %>% pull(host_genus)


#### Make genera found in bacterial analysis bold

## Get IBS-assoicated genera 
top_genus <- gsub('_.*', '', top_bacteria$X)

## Add bold labels to genus name if this genus is also in IBS-associated bacteria
host_OR_counts$host_genus_bold <- ifelse(host_OR_counts$host_genus %in% top_genus, paste0("**", host_OR_counts$host_genus, "**"), host_OR_counts$host_genus)

## Enforce ordering after updating text labels
host_OR_counts$host_plot <- factor(
  host_OR_counts$host_genus_bold,
  levels = unique(host_OR_counts$host_genus_bold[order(match(host_OR_counts$host_genus, ordering))])
)

plot <- ggplot(host_OR_counts, aes(x = host_plot, y = n, fill = OR_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "steelblue")) +
  labs(x = "Predicted Bacterial Host", y = "Count", fill = "vOTU ssociation\n with IBS") +
  ggtitle("Host Predictions \nof IBS-Associated vOTUs") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.y = ggtext::element_markdown(size = 25, margin = margin(r = 10)),
    axis.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 31),
    plot.title = element_text(hjust = 0.5, size = 32),
    legend.title = element_text(size = 21),
    legend.text = element_text(size = 20),
    plot.margin = unit(c(2, 1, 0.1, 1), "lines")
  )




# Host-adjusted  ----------------------------------------------------------

## Load host-adjusted results
host_adj <- read.csv("data_deposition/gut_paper/gut_paper_results/phage_results/phage_lmm_host_adjusted.csv")

## Change rownames to virus name for merging
host_adj <- host_adj %>% rename(Virus = X)

## Add analysis-type for plotting
host_adj$Analysis <- 'Host-adjusted'


## Add host genus labels to host_adj
host_labs <- c()
for (i in host_adj$Virus) {
  lab <- paste0(gsub('.*g__', '', unique(hosts_ibs$host_genus[hosts_ibs$sudo_phage_ID == i])), collapse =', ')
  host_labs <- c(host_labs, lab)
}
host_adj$Bacteria <- host_labs


## How many associations no longer significant after adjusting for host?
sum(host_adj$PVAL > 0.05) # 0

## Keep phages from main analysis that were also in host-adjusted analysis
host_main <- res_df[res_df$Virus %in% host_adj$Virus, ]
host_main$Analysis <- 'Overall'
host_main <- merge(host_main, host_adj[c('Virus', 'Bacteria')], by = 'Virus')

cols <- c('Virus', 'Bacteria', 'OR', 'Lower', 'Upper', "PVAL", 'Analysis')
host_adj_df <- rbind(host_main[cols], host_adj[cols])

# Update labels 
host_adj_df$label <- paste0('(', host_adj_df$Bacteria, ') ', host_adj_df$Virus)
host_adj_df$label <- gsub('_', ' ', host_adj_df$label)

# Sort by OR original
order <- host_adj_df[host_adj_df$Analysis == 'Overall', ] %>% arrange(OR) %>% pull(label)
host_adj_df$label <- factor(host_adj_df$label, levels = order) 
host_adj_df <- host_adj_df %>% arrange(OR)

# Format p-value for plot
host_adj_df$PVAL <- ifelse(host_adj_df$PVAL < 0.01, formatC(host_adj_df$PVAL, format = "e", digits = 2), signif(host_adj_df$PVAL, 2))


position <- position_dodge2(width =  0.6)  
plot <- ggplot(host_adj_df, aes(y = label, x = OR, color = Analysis)) +
  geom_point(size = 3, position = position) + 
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height =  0.6, position = position) +
  geom_vline(xintercept = 1, color = "#0073C2FF", linetype = "dashed", size = 1, alpha = 0.5) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.4))  # 20% padding to the right
  ) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Phage (predicted host)") + 
  ggtitle('IBS-associated vOTUs With and Without Adjustment\n for Predicted Host Genus Relative Abundance\n') + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.text.y = element_text(size = 29, colour = "black", face = "italic"),
        axis.title.x = element_text(size = 29, colour = "black"), 
        axis.title.y = element_text(size = 33, colour = "black"), 
        plot.title = element_text(hjust = 0.92, size = 42),
        plot.margin = unit(c(1, 11, 0.1, 3), "lines"),
        strip.text.y = element_text(size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 31),
        legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 1)) +
  coord_cartesian(clip = "off") + 
  geom_text(data = host_adj_df %>% filter(Analysis == "Host-adjusted"),
            aes(label = paste0('p=', PVAL)),
            x = (max(host_adj_df$Upper)+0.06), 
            hjust = 0, 
            size = 11, 
            color = "black")


# AMG analysis -------------------------------------------------

amg <- read.csv("data_deposition/gut_paper/gut_paper_datasets/metagenome_data/phage_data/Database_independent_approach_revised/viral contigs metadata/all_amg_summary_final_contigs.tsv", sep = '\t')

## Exclude low confidence AMG assignments 
amg <- amg[amg$auxiliary_score == 1, ] # "1 represents a gene that is confidently virally encoded and a score of 4 or 5 represents a gene that users should take caution in treating as a viral gene" (https://pmc.ncbi.nlm.nih.gov/articles/PMC7498326/)

## Keep only those in our dataset & add pseudo phage ID
amg_sig <- merge(amg, unique(hosts_ibs[, c('Virus', 'sudo_phage_ID')]), by.x = 'gene', by.y = 'Virus')

## Keep only those labelled as potential AMGs
amg_sig <- amg_sig[amg_sig$potential_amg == 'True', ]

## How many AMGs are present among IBS-associated VCs?
length(unique(amg_sig$gene_id))                         # 2 unique gene IDs
length(unique(amg_sig$module))                          # 2 modules

## From how many unique phages do these AMGs come?
length(unique(amg_sig$gene))                            # Two unique phages

## What was direction of change of relative abundance phage in results dataset?
res_df$OR[res_df$Virus %in% amg_sig$sudo_phage_ID]      # Both were lower in IBS

## What are the names of the represented AMGs?
table(amg_sig$gene_id)
# GH13    M26 
# 1       1       

## What are the descriptions of these genes?
unique(amg_sig$gene_description)
# [1] "Catlytic type: Metallo; cleaves the heavy chain of human IgA1 at the Pro227-Thr228 bond; tightly associated with the bacterial cell surface"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
# [2] "GH13 alpha-amylase (EC 3.2.1.1); pullulanase (EC 3.2.1.41); cyclomaltodextrin glucanotransferase (EC 2.4.1.19); cyclomaltodextrinase (EC 3.2.1.54); trehalose-6-phosphate hydrolase (EC 3.2.1.93); oligo-alpha-glucosidase (EC 3.2.1.10); maltogenic amylase (EC 3.2.1.133); neopullulanase (EC 3.2.1.135); alpha-glucosidase (EC 3.2.1.20); maltotetraose-forming alpha-amylase (EC 3.2.1.60); isoamylase (EC 3.2.1.68); glucodextranase (EC 3.2.1.70); maltohexaose-forming alpha-amylase (EC 3.2.1.98); maltotriose-forming alpha-amylase (EC 3.2.1.116); branching enzyme (EC 2.4.1.18); trehalose synthase (EC 5.4.99.16); 4-alpha-glucanotransferase (EC 2.4.1.25); maltopentaose-forming alpha-amylase (EC 3.2.1.-) ; amylosucrase (EC 2.4.1.4) ; sucrose phosphorylase (EC 2.4.1.7); malto-oligosyltrehalose trehalohydrolase (EC 3.2.1.141); isomaltulose synthase (EC 5.4.99.11); malto-oligosyltrehalose synthase (EC 5.4.99.15); amylo-alpha-1,6-glucosidase (EC 3.2.1.33); alpha-1,4-glucan: phosphate alpha-maltosyltransferase (EC 2.4.99.16); 6'-P-sucrose phosphorylase (EC 2.4.1.-); amino acid transporter"  


## Which modules are present?
table(amg_sig$module)
# Endopeptidases       Glycoside Hydrolases 
# 1                    1 

## Which values are present in "sheets"?
table(amg_sig$sheet)
# Unknown               carbon utilization   Organic Nitrogen 
# 12                    1                    1 

table(amg_sig$header)
# Unknown               CAZY                 Peptidase 
# 12                    1                    1 

table(amg_sig$subheader)
# Unknown               Starch Backbone Cleavage 
# 13                    1 

# Only two of these - those which are not "Unknown" in the previous 4 tables - are "True" for "potential_amg"

red <- amg_sig[amg_sig$module != '', c('gene', "gene_description", 'module', 'sheet', 'header', 'subheader')]













