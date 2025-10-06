# laguna-lake-ARG-Analysis
Laguna Lake, the Philippines’ largest freshwater lake, harbors antibiotic-resistant bacteria, yet data on antibiotic resistance genes (ARGs) remain scarce. This repository serves as documentation of the data analysis of a metagenomic study investigating the diversity, prevalence, and mobility of ARGs in the lake’s West Bay.

## documentation files
- `ARG-Quantification-Workflow.md` : Quantifies ARGs at type and subtype level, and generates a heatmap and a stacked bar chart.
- `ARC-Profiling-Workflow.md` : Analyzes mobility of ARGs based on the mobility characteristics (plasmid and MGE association) of the contigs that contain them.
## data availability
Currently, the data associated with this analysis is not yet uploaded to any public repository.
## contents
- `ARC_faa/` and `ARC_fasta/` : contain fasta files of antibiotic-resistant contigs (ARCs)
- `BWA_output/` : contains the resulting bam files of BWA-aligned short reads to contigs
- `Diamond_SARG_output/` : DIAMOND alignments of contigs to SARG database
- `Diamond_mobileOG_output/` : DIAMOND alignments of contigs to MobileOG database
- `Plasflow/` : Plasflow classifications of contigs
- `args-oap_output/` : contains output files of ARGs-OAP run
- `clean_reads/` : contains the clean reverse and forward reads
- `contig_coverage_stats/` : contains the coverage statistics of contigs
- `contigs/` : contains assembled contigs
- `contigs_500bp/` : contains contigs >= 500bp
- `databases/` : contains the SARG and MobileOG databases used for the analysis
- `figures/` : contains output plots
- `predicted_ORFs/` : contains the open reading frames predicted by Prodigal
- `raw_reads/` : contains the raw reverse and forward reads
- `tables/` : contains the output tables - supplementary data
