get_multiplex_annotations.sh - master script

#fetches all mutations from central database
#creates vcf file
#if have genomic locations - get all COSMIC annotations, or else corrects the amino acid change
#if no genomic locations - fetches possible genomic location as vcf file
#maps annotations back to mutations csv file
