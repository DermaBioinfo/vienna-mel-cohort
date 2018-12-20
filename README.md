# vienna-mel-cohort

A simple R shiny App summarizing the basic characteristics of the Vienna Melanoma cohort based on the published supplementary data.

## Specification

  * Displays one table with all melanoma samples
  * One column for every tested mutated gene (for example "NF1", "BRAF" etc)
    * First version only contains BRAF, NRAS and NF1
    * Should have option to easily add additional genes
  * User interace with one drop-down per gene with 3 values ("ignore", "mutated", "non-mutated")
