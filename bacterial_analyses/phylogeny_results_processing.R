# Explanation -------------------------------------------------------------

# Uses the results generated from phylogeny_analysis.py to identify the 
# nearest known relative of unknown bacterial species associated with IBS



# Load tree and results dataframes ----------------------------------------

csv_files <- list.files('data_deposition/gut_paper/gut_paper_results/bacterial_results/phylogeny_tree', pattern = "//.csv$", full.names = TRUE)

# Read each file and assign it to an object named after the file
dfs <- setNames(
  lapply(csv_files, read.csv),
  tools::file_path_sans_ext(basename(csv_files))
)


# Get closest known relative ----------------------------------------------

# If all top relatives are unknown, return name for closest inspection 

top_known <- c()
for (i in names(dfs)) {
  df <- dfs[[i]][-1]
  
  # Account for multiple species entries per SGB
  df$full_name <- sub(',.*', '', df$full_name)
  
  # Tidy up species name and add genus
  df$species <- sub(".*s__", "", df$full_name)
  df$genus <- sub(".*g__([^|]+).*", "\\1", df$full_name)
  df$family <- sub(".*f__([^|]+).*", "\\1", df$full_name)
  df <- df[, c('SGB', 'full_name', 'species', 'genus', 'family',  'distance')]
  
  cat('\n', i)
  print(df)
  
  # Keep all remaining rows
  top_known[[i]] <- df
  
  # Save as object
  assign(i, df, envir = .GlobalEnv)
}







