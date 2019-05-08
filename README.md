# VCFAnnotation
Simple VCF annotation tool Implemented as RShiny app

Based on *Ensembl* Release 75 and ExAC database

Input: a VCF file

Each variant is annotated with the following pieces of information:

1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API (API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional optional information from ExAC.

An example vcf file data/Challenge_data (1).vcf is provided for testing purpose. The annotated vcf downloaded from the tool is also provided at data/annotated_vcf.vcf. It will take ~30 minutes to run the whole file (~7000 variants). 

To run from zip: unzip and open the "app.R" in RStudio. Click "Run App" button to run. 

To run from Github: https://github.com/hoyden0329/VCFAnnotation.git

`runGitHub( "VCFAnnotation", "hoyden0329")` 

To run from Shinyapps.io (free account, 25 hours usage limit per month; memory limit, cannot run large vcf files):  https://xiuhuang.shinyapps.io/VCFAnnotation/