# Load necessary packages
library(tidyr)
library(dplyr)
library(ggplot2)

#### Profiling by subtype; normalized to cell count
# Import args-oap output file: metadata.txt
args_oap_metadata <- read.table('args-oap_output/metadata.txt', header=TRUE, row.names=1)

# Import args-oap output file: normalized_cell.subtype
normalized_cell.subtype <- as.data.frame(read.delim('args-oap_output/normalized_cell.subtype.txt', header=TRUE))
colnames(normalized_cell.subtype) <- sub("_.*", "", colnames(normalized_cell.subtype)) # Clean column names

# Separate type and subtype names
ARG_subtype_df <- normalized_cell.subtype %>%
  separate(col = subtype, into = c("Type", "Subtype"), sep = "__") %>%
  select(Type, Subtype, everything())

# Add OVERALL abundance column. OVERALL represents the cumulative abundance of ARG types from all samples normalized to the total cell copy number
ARG_subtype_df$OVERALL <- (ARG_subtype_df$AUGLS01*args_oap_metadata[[1, "nCell"]] +
                             ARG_subtype_df$AUGLS05*args_oap_metadata[[2, "nCell"]] +
                             ARG_subtype_df$AUGLS16*args_oap_metadata[[3, "nCell"]] +
                             ARG_subtype_df$SEPLS01*args_oap_metadata[[4, "nCell"]] +
                             ARG_subtype_df$SEPLS05*args_oap_metadata[[5, "nCell"]] +
                             ARG_subtype_df$SEPLS16*args_oap_metadata[[6, "nCell"]]
                            ) / sum(args_oap_metadata[["nCell"]])

# List of ARG types (overall abundance > 0.009 cpc). Check ARG_type_heatmap.R to see abundant ARG types.
Abundant_ARGs <- c("beta_lactam", "bacitracin", "multidrug", "chloramphenicol", "polymyxin", "tetracycline")
# Filter abundant ARG types
Abundant_ARGs_df <- ARG_subtype_df %>%
  filter(grepl(paste(Abundant_ARGs, collapse = "|"), Type)) %>%
  # Set order of ARG types
  mutate(Type = factor(Type, levels = c(
    "beta_lactam", "bacitracin", "multidrug", "chloramphenicol", "polymyxin", "tetracycline"
  )))


# Group ARGs based on subtypes and compute initial total abundance per subtype
Abundant_ARGs_grouped <- Abundant_ARGs_df %>%
  group_by(Type, Subtype) %>%
  summarise(Total_Abundance = sum(OVERALL), .groups = "drop") %>%
  group_by(Type) %>%
  mutate(
    Relative_Abundance = Total_Abundance / sum(Total_Abundance),
    subtype_substr = substr(Subtype, 1, 3)
    ) %>%
  # manually group subtypes based on ARG families in CARD (https://card.mcmaster.ca/)
  mutate(subtype_group = case_when(
    # for bacitracin resistance, keep subtypes as is since they are all from different gene families
    Type == "bacitracin" ~ Subtype,
    # for beta lactam resistance, the default grouping is the first 3 letters of the subtype (dropping the number)
    Type == "beta_lactam" & Subtype == "Escherichia coli ampC" ~ "ampC", # ampC-type beta-lactamase
    Type == "beta_lactam" & Subtype == "Klebsiella pneumoniae OmpK37" ~ "OmpK37",
    Type == "beta_lactam" & Subtype == "Other class A beta-lactamase" ~ Subtype,
    Type == "beta_lactam" & Subtype == "Other class C beta-lactamase" ~ Subtype,
    Type == "beta_lactam" ~ subtype_substr, # default for beta lactams
    # For multidrug resistance, subtypes are grouped based on major gene family
    Type == "multidrug" & Subtype == "Enterobacter cloacae acrA" ~ "RND type drug efflux",
    Type == "multidrug" & Subtype == "Escherichia coli acrA" ~ "RND type drug efflux",
    Type == "multidrug" & Subtype == "Escherichia coli emrE" ~ "SMR type drug efflux",
    Type == "multidrug" & Subtype == "Escherichia coli mdfA" ~ "MFS type drug efflux",
    Type == "multidrug" & subtype_substr == "Mex" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "Mux" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "Opm" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "Opr" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "Ran" ~ "ABC type drug efflux",
    Type == "multidrug" & Subtype == "abeM" ~ "MATE transporter",
    Type == "multidrug" & subtype_substr == "acr" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "ade" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "amr" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "ceo" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "efp" ~ "MFS type drug efflux",
    Type == "multidrug" & subtype_substr == "efr" ~ "ABC type drug efflux",
    Type == "multidrug" & subtype_substr == "eme" ~ "MATE transporter",
    Type == "multidrug" & subtype_substr == "emr" ~ "MFS type drug efflux",
    Type == "multidrug" & subtype_substr == "mds" ~ "RND type drug efflux",
    Type == "multidrug" & Subtype == "mdtE" ~ "RND type drug efflux",
    Type == "multidrug" & Subtype == "mdtF" ~ "RND type drug efflux",
    Type == "multidrug" & Subtype == "mdtH" ~ "MFS type drug efflux",
    Type == "multidrug" & Subtype == "mdtK" ~ "MATE transporter",
    Type == "multidrug" & subtype_substr == "mdt" ~ "MFS type drug efflux",
    Type == "multidrug" & subtype_substr == "mex" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "msb" ~ "ABC type drug efflux",
    Type == "multidrug" & subtype_substr == "mtr" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "oqx" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "pmr" ~ "MFS type drug efflux",
    Type == "multidrug" & Subtype == "qacE" ~ "MFS type drug efflux",
    Type == "multidrug" & Subtype == "qacEdelta1" ~ "MFS type drug efflux",
    Type == "multidrug" & Subtype == "qacH" ~ "SMR type drug efflux",
    Type == "multidrug" & Subtype == "sdeY" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "sme" ~ "RND type drug efflux",
    Type == "multidrug" & subtype_substr == "tap" ~ "MFS type drug efflux",
    # for chloramphenicol resistance, grouped based on gene family except for catI (~93% rel abundance)
    Type == "chloramphenicol" & Subtype == "catI" ~ "catI",
    Type == "chloramphenicol" & Subtype == "plasmid-encoded cat (pp-cat)" ~ "Other cat",
    Type == "chloramphenicol" & subtype_substr == "cat" ~ "Other cat", # default for cat except catI 
    Type == "chloramphenicol" & Subtype == "catI" ~ Subtype,
    Type == "chloramphenicol" & Subtype == "Salmonella enterica cmlA" ~ "cml",
    Type == "chloramphenicol" & subtype_substr == "cml" ~ subtype_substr,
    # for polymyxin resistance, most are pmr phosphoethanolamine transferase family.
    # but since ugd accounts for (~71%), it will be kept as is except for other pmr genes
    Type == "polymyxin" & Subtype == "arnA" ~ "Other pmr",
    Type == "polymyxin" & Subtype == "eptA" ~ "Other pmr",
    Type == "polymyxin" & Subtype == "pmrF" ~ "Other pmr",
    Type == "polymyxin" & subtype_substr == "mcr" ~ subtype_substr, # for mcr-5.1 and -5.2
    # for rosA and rosB, they will be considered as rosAB as they make up the two-component system
    Type == "polymyxin" & Subtype == "rosA" ~ "rosAB",
    Type == "polymyxin" & Subtype == "rosB" ~ "rosAB",
    # for tetracycline resistance, the subtypes are kept as is
    TRUE ~ Subtype # default for all others that are not renamed
  ))

# Compute abundance by subtype group
ARG_group_abundance <- Abundant_ARGs_grouped %>%
  group_by(Type, subtype_group) %>%
  summarise(Group_Abundance = sum(Total_Abundance), .groups = "drop") %>%
  group_by(Type) %>%
  mutate(Relative_Abundance = Group_Abundance / sum(Group_Abundance)) %>%
  # recategorize groups with < 5% abundance as Other 
  mutate(subtype_group = ifelse(Relative_Abundance < 0.05, "Other", subtype_group)) %>%
  group_by(Type, subtype_group) %>%
  summarise(Group_Abundance = sum(Group_Abundance), .groups = "drop") %>%
  group_by(Type) %>%
  mutate(Relative_Abundance = Group_Abundance / sum(Group_Abundance)) %>%
  # Reorder by ARG type and groups for stacked bar chart (Other on top, most abundant at the bottom)
  mutate(subtype_group = factor(subtype_group, levels = c(
    "Other",
    "TEM",
    "bcrA", "bacA",
    "RND type drug efflux", "ABC type drug efflux", "MFS type drug efflux",
    "Other cat", "catI",
    "rosAB", "Other pmr","ugd",
    "tet(B)", "tet(C)"
  )))

## plot
# set colors
group_colors <- c(
  "Other" = "azure2",
  "TEM" = "brown1",
  "bacA" = "deepskyblue4", "bcrA" = "deepskyblue3",
  "MFS type drug efflux" = "darkolivegreen4", "ABC type drug efflux" = "darkolivegreen3", "RND type drug efflux" = "darkolivegreen1",
  "catI"                 = "coral3", "Other cat" = "coral",
  "ugd" = "darkgoldenrod4", "Other pmr" = "darkgoldenrod3", "rosAB" = "darkgoldenrod1",
  "tet(C)" = "darkseagreen3", "tet(B)" = "darkseagreen1"
)

# create stacked bar chart
ggplot(ARG_group_abundance, aes(x = Type, y = Relative_Abundance, fill = subtype_group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "ARG Type",
    y = "Relative Abundance"
  ) +
  guides(fill = guide_legend(title = NULL)) +
  theme_minimal(base_family = "sans", base_size = 12) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0),
        panel.background = element_blank())+
  scale_fill_manual(values = group_colors)


