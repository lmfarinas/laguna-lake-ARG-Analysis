# Load necessary packages
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(scales)

# List of sample names
samples <- c("AUGLS01", "AUGLS05", "AUGLS16", "SEPLS01", "SEPLS05", "SEPLS16")

# Read pertinent files
SargDB <- read.delim("databases/SARG_v3.2.1-S/4.SARG_v3.2_20220917_Short_subdatabase_structure.txt", header = TRUE, stringsAsFactors = FALSE) # sarg database structure
ARGs_metadata <- read.delim("args-oap_output/metadata.txt", header = TRUE, stringsAsFactors = FALSE) # metadata from args-oap analysis
mobileOG1.6 <- read.delim(file = "databases/mobileOG/mobileOG-db-beatrix-1.6-All.csv", sep = ",") # mobileOG1.6 database structure
plasflow_prob <- read.delim("Plasflow/ARCs_plasflow_prob_table.tabular", header = TRUE, stringsAsFactors = FALSE) # plasflow probability table

#### Calculate abundance of ARCs
# extract data for ncell data
ncell_data <- ARGs_metadata[,-c(2,3)]
ncell_data$Sample <- sub("_.*", "", ARGs_metadata$Sample)
# Calculate overall ncell
overall_ncell <- sum(ncell_data$nCell)

# Function to read coverage stats and calculate ARC abundance for a sample file
getCoverageStats <- function(Sample) {
  # Read coverage stats
  cov_stats <- read.delim(paste0("contig_coverage_stats/", Sample, "_500bpContigs_cov_stats.txt"), header = TRUE, stringsAsFactors = FALSE)
  
  # Filter for ARCs
  ARC_hits <- read.delim(paste0("Diamond_SARG_output/", Sample, "_dmnd_sarg.txt"), header = FALSE, stringsAsFactors = FALSE) # read diamond hits with SARG db
  ARC_nodes <- unique(sub("_[^_]*$", "", ARC_hits[[1]]))        # get the unique ARCs
  ARC_stats <- cov_stats[cov_stats$X.rname %in% ARC_nodes, ]    # Only get stats of ARCs
  
  ARC_stats <- ARC_stats %>%
    rename(Contig_ID = X.rname) %>%   # Rename the column
    mutate(Sample = Sample) %>%       # Add the Sample column
    select(Sample, everything()) %>%  # Reorder to make Sample the first column
    mutate(Contig_ID = sub("^(([^_]*_[^_]*)).*", "\\1", Contig_ID)) %>%    # drop length and coverage in Contig_ID
    mutate(Contig_ID = paste0(Sample, "_", Contig_ID))                     # include sample name to Contig_ID
  
  # Calculate abundance
  ncell <- ncell_data[ncell_data$Sample==Sample, 2]                               # Get ncell for the sample
  ARC_stats$Abundance_cpc <- ARC_stats$numreads * 150 / ARC_stats$endpos / ncell  # formula for abundance of an ARC
  
  # clear row names
  rownames(ARC_stats) <- NULL
  
  return(ARC_stats)
}

# Apply the function to each sample and combine results
all_cov_stats <- do.call(rbind, lapply(samples, getCoverageStats)) # contains contig coverage statistics

#### Classify ARCs based on plasflow probability table
plasflow_classification <- plasflow_prob %>%
  mutate(
    Contig_ID = str_replace(contig_name, "^((?:[^_]+_){2}[^_]+).*", "\\1"),
    PlasFlow_classification = str_extract(label, "^[^\\.]+")
  ) %>%
  select(Contig_ID, PlasFlow_classification)

#### Annotate ARCs with their top ARG and MGE hits
# Function to Annotate ARCs with ARG and MGE hits
getARGMGEHits <- function(Sample) {
  
  # Read hits
  ARG_hits <- read.delim(paste0("Diamond_SARG_output/", Sample, "_dmnd_sarg.txt"), header = FALSE, stringsAsFactors = FALSE)              # read diamond hits with SARG db
  MGE_hits <- read.delim(paste0("Diamond_mobileOG_output/", Sample, "_ARCs_dmnd_mobileOG.txt"), header = FALSE, stringsAsFactors = FALSE) # read diamond hits with MobileOG db
  
  # Add column names
  colnames(ARG_hits) <- c("qseqid", "sseqid", "salltitles", "pident", "length", "mismatch", 
                      "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  colnames(MGE_hits) <- c("qseqid", "sseqid", "salltitles", "pident", "length", "mismatch", 
                          "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  
  ARG_hits <- ARG_hits %>%
    mutate(ORF = sub(".*_(\\d+)$", "\\1", qseqid)) %>%            # Add the ORF column
    mutate(database = "SARG") %>%                                 # Add the database column
    mutate(qseqid = sub("^(([^_]*_[^_]*)).*", "\\1", qseqid)) %>% # drop length and coverage in qseqid
    mutate(qseqid = paste0(Sample, "_", qseqid)) %>%              # include sample name to qseqid
    rename(Contig_ID = qseqid) %>%                                # Rename the column
    mutate(Sample = Sample) %>%                                   # Add the Sample column
    select(-salltitles) %>%                                       # Remove salltitles column
    select(Sample, Contig_ID, ORF, database, everything())        # Reorder columns
  
  MGE_hits <- MGE_hits %>%
    mutate(ORF = sub(".*_(\\d+)$", "\\1", qseqid)) %>%            # Add the ORF column
    mutate(database = "MobileOG") %>%                             # Add the database column
    mutate(qseqid = sub("^(([^_]*_[^_]*)).*", "\\1", qseqid)) %>% # drop length and coverage in qseqid
    mutate(qseqid = paste0(Sample, "_", qseqid)) %>%              # include sample name to qseqid
    rename(Contig_ID = qseqid) %>%                                # Rename the column
    mutate(Sample = Sample) %>%                                   # Add the Sample column
    select(-salltitles) %>%                                       # Remove salltitles column
    select(Sample, Contig_ID, ORF, database, everything())        # Reorder columns
  
  # Merge dataframes
  ARGMGE_hits <- rbind(ARG_hits, MGE_hits)
  
  # Sort rows
  ARGMGE_hits <- ARGMGE_hits %>%
    mutate(Contig_num = as.numeric(str_extract(Contig_ID, "(?<=NODE_)\\d+"))) %>%   # Add a sorting column
    arrange(Contig_num, ORF) %>%
    select(-Contig_num)   # remove sorting column

  # clear row names
  rownames(ARGMGE_hits) <- NULL
  
  return(ARGMGE_hits)
}

# Apply the function to each sample and combine results
ARC_ARGMGE_hits <- do.call(rbind, lapply(samples, getARGMGEHits)) # contains top Diamond hits (with SargDB and MobileOG) per ORF of ARCs

# Annotate hits per ORF based on sseqid
Annotated_ARCs_hits <- ARC_ARGMGE_hits %>% # contains annotation from the databases (ARG type and subtype, and MGE category)
  # Add ARG type and subtype using SargDB structure
  left_join(SargDB %>% select(`SARG.Seq.ID`, Type, Subtype),
            by = c("sseqid" = "SARG.Seq.ID")) %>%
  rename(ARG_type = Type, ARG_subtype = Subtype) %>%
  mutate(ARG_subtype = sapply(strsplit(ARG_subtype, "__"), `[`, 2)) %>%
  # Add MGE category using sseqid string
  mutate(MGE_category = if_else(
      database == "MobileOG",sapply(strsplit(sseqid, "\\|"), `[`, 4), NA_character_)) %>%
  # Use the mobileOG1.6 db structure file to annotate sseqids with no info
  mutate(sseqid_1 = if_else(
    database == "MobileOG" & is.na(MGE_category), sapply(strsplit(sseqid, "\\|"), `[`, 1), NA_character_)) %>%
  left_join(
    mobileOG1.6 %>% select(mobileOG.Entry.Name, Major.mobileOG.Category),
    by = c("sseqid_1" = "mobileOG.Entry.Name")
  ) %>%
  mutate(MGE_category = if_else(
    !is.na(Major.mobileOG.Category), Major.mobileOG.Category, MGE_category)) %>%
  select(-sseqid_1, -Major.mobileOG.Category)   # remove sseqid_1 and Major.mobileOG.Category columns

#### Abbreviate ARG types and MGE categories
# Define a named vector for ARG type abbreviations
arg_abbr <- c(
  tetracycline = "tet",
  beta_lactam = "bla",
  chloramphenicol = "chl",
  sulfonamide = "sul",
  multidrug = "mul",
  novobiocin = "nov",
  polymyxin = "pol",
  bacitracin = "bac",
  trimethoprim = "tri",
  `macrolide-lincosamide-streptogramin` = "mls",
  pleuromutilin_tiamulin = "plt",
  mupirocin = "mup",
  aminoglycoside = "ami",
  other_peptide_antibiotics = "opa"
)
# Define a named vector for MGE category abbreviations
mge_abbr <- c(
  "replication/recombination/repair" = "RRR",
  "integration/excision" = "IE",
  transfer = "T",
  phage = "P",
  "stability/transfer/defense" = "STD"
)

# get selected columns and abbreviate ARG types
ARCs_arg_abbr <- Annotated_ARCs_hits %>%
  filter(!is.na(ARG_type)) %>% # keep only the rows with ARG hits
  select(Contig_ID, ARG_type) %>% # Keep only selected columns
  mutate(ARG_type = recode(ARG_type, !!!arg_abbr)) %>% # Replace ARG_type with abbreviations
  # Collapse distinct ARG types per Contig_ID
  group_by(Contig_ID) %>%
  summarize(
    ARG_type = paste(unique(ARG_type), collapse = "-"),
    .groups = "drop"
  )

# get selected columns and abbreviate MGE categories
ARCs_mge_abbr <- Annotated_ARCs_hits %>%
  filter(!is.na(MGE_category)) %>% # keep only the rows with MGE hits
  select(Contig_ID, MGE_category) %>% # Keep only selected columns
  mutate(MGE_category = recode(MGE_category, !!!mge_abbr)) %>% # Replace MGE_category with abbreviations
  # Collapse distinct ARG types per Contig_ID
  group_by(Contig_ID) %>%
  summarize(
    MGE_category = paste(unique(MGE_category), collapse = "-"),
    .groups = "drop"
  )

#### Summarize ARCs annotation
# create table for abundance, ARG, MGE, plasflow classification
ARCs_summary <- all_cov_stats %>%
  select(Sample, Contig_ID, Abundance_cpc) %>% # keep pertinent columns
  left_join(ARCs_arg_abbr, by = "Contig_ID") %>%  # Add ARG type
  left_join(ARCs_mge_abbr, by = "Contig_ID") %>%  # Add MGE category
  left_join(plasflow_classification, by = "Contig_ID")  # Add PlasFlow classification

# List of ARCs with associated ARGs > 0.001 cpc based on ARGs-OAP quantification
arc_0.001 <- c("bla", "bac", "mul", "chl", "pol", "tet", "mup", "mls", "sul")

ARCs_summary_0.001 <- ARCs_summary %>%
  separate_rows(ARG_type, sep = "-") %>%    # Split multiple ARGs into rows
  filter(ARG_type %in% arc_0.001) %>%       # Keep only target ARGs using arc_0.001
  left_join(ncell_data, by = "Sample") %>%  # Add ncell values
  mutate(UN_cpc = Abundance_cpc*nCell)      # Compute unnormalized abundances
  
### Summarize plasmid association of ARCs
plasmid_summary <- ARCs_summary_0.001 %>%
  mutate(Plasmid_Status = ifelse(PlasFlow_classification == "plasmid", "Plasmid", "Non-Plasmid")) %>%   # Classify either as Plasmid or Non-Plasmid
  group_by(ARG_type, Plasmid_Status) %>%                          # Group rows by their ARG_type and Plasmid_Status
  summarise(Total_Abundance = sum(UN_cpc), .groups = "drop") %>%  # Calculate total unnormalized abundance
  group_by(ARG_type) %>%                                          # Group by ARG_type
  mutate(Relative_Abundance = Total_Abundance / sum(Total_Abundance)) # calculate relative abundance

# Calculate overall ratio of Plasmid vs Non-plasmid (not just for the selected ARGs)
plasmid_overall <- ARCs_summary  %>%
  left_join(ncell_data, by = "Sample") %>%  # Add ncell values
  mutate(Plasmid_Status = ifelse(PlasFlow_classification == "plasmid", "Plasmid", "Non-Plasmid")) %>%   # Classify either as Plasmid or Non-Plasmid
  mutate(UN_cpc = Abundance_cpc*nCell) %>%  # Compute unnormalized abundances
  group_by(Plasmid_Status) %>%              # Group by plasmid status
  summarise(Total_Abundance = sum(UN_cpc), .groups = "drop") %>%    # Calculate total unnormalized abundance
  mutate(ARG_type = "ALL") %>%              # Change ARG_type to "ALL"
  group_by(ARG_type) %>%              
  mutate(Relative_Abundance = Total_Abundance / sum(Total_Abundance)) # Calculate relative abundance

# Combine plasmid association by ARG type with overall
plasmid_combined <- bind_rows(plasmid_summary, plasmid_overall)

# Set order of factors
plasmid_combined$ARG_type <- factor(plasmid_combined$ARG_type,
                                               levels = c("ALL", "bla", "bac", "mul", "chl", "pol", "tet", "mup", "mls", "sul"))

# Plot
ggplot(plasmid_combined, aes(x = ARG_type, y = Relative_Abundance, fill = Plasmid_Status)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "ARC",
    y = "Relative Abundance"
  ) +
  geom_text(
    data = plasmid_combined %>% filter(Plasmid_Status == "noPlasmid"),
    aes(label = percent(Relative_Abundance, accuracy = 1)),
    position = position_stack(vjust = 0.5), # center of Plasmid segment
    color = "white",
    size = 3.5,
    family = "sans"
  ) +
  guides(fill = guide_legend(title = NULL)) +
  theme_minimal(base_family = "sans", base_size = 15) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values = c("Plasmid" = "#8B0000", "Non-Plasmid" = "grey"))

### Summarize MGE association of ARCs
MGE_summary <- ARCs_summary_0.001 %>%
  mutate(MGE_category = ifelse(is.na(MGE_category), "no MGE gene", MGE_category)) %>%    # annotate ARCs with no MGE as no MGE gene
  group_by(ARG_type, MGE_category) %>%                          # Group rows by their ARG_type and MGE_category
  mutate(MGE_category = if_else(MGE_category == "IE-RRR", "RRR-IE", MGE_category)) %>%    # Reorder MGEs to match with others
  summarise(Total_Abundance = sum(UN_cpc), .groups = "drop") %>%  # Calculate total unnormalized abundance
  group_by(ARG_type) %>%                                          # Group by ARG_type
  mutate(Relative_Abundance = Total_Abundance / sum(Total_Abundance)) # calculate relative abundance

# Calculate overall ratio of Plasmid vs Non-plasmid (not just for the selected ARGs)
MGE_overall <- ARCs_summary  %>%
  left_join(ncell_data, by = "Sample") %>%  # Add ncell values
  mutate(MGE_category = ifelse(is.na(MGE_category), "no MGE gene", "w/ MGE gene")) %>%    # annotate ARCs with no MGE as no MGE gene
  mutate(UN_cpc = Abundance_cpc*nCell) %>%  # Compute unnormalized abundances
  mutate(ARG_type = "ALL") %>%              # Change ARG_type to "ALL"
  group_by(ARG_type, MGE_category) %>%
  summarise(Total_Abundance = sum(UN_cpc), .groups = "drop") %>%      # Compute total unnormalized abundances
  mutate(Relative_Abundance = Total_Abundance / sum(Total_Abundance)) # Compute relative abundance

# Combine MGE association by ARG type with overall
MGE_combined <- bind_rows(MGE_summary, MGE_overall)

# Set order of factors
MGE_combined$ARG_type <- factor(MGE_combined$ARG_type,
                                    levels = c("ALL", "bla", "bac", "mul", "chl", "pol", "tet", "mup", "mls", "sul"))

# Plot
MGE_colors <- c(
  "no MGE gene" = "grey",
  "w/ MGE gene" = "brown1",
  "RRR"         = "deepskyblue4",
  "RRR-IE"          = "darkolivegreen3",
  "T"          = "coral",
  "IE"        = "darkgoldenrod4",
  "P"          = "darkgoldenrod1"
)

MGE_combined <- MGE_combined %>%
  mutate(MGE_category = factor(MGE_category, levels = c(
    "no MGE gene", "w/ MGE gene", "RRR-IE", "RRR", "IE", "T", "P"))) # Set order of factors

ggplot(MGE_combined, aes(x = ARG_type, y = Relative_Abundance, fill = MGE_category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "ARC",
    y = "Relative Abundance"
  ) +
  geom_text(
    data = MGE_combined %>% filter(MGE_category == "w/in MGE gene"),
    aes(label = percent(Relative_Abundance, accuracy = 1)),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3.5,
    family = "sans"
  ) +
  guides(fill = guide_legend(title = NULL)) +
  theme_minimal(base_family = "sans", base_size = 15) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values = MGE_colors)
