##calculate candidate gene of Fst
bedtools intersect -a fst_candidate_region.txt -b B_ZM4_bed_region.txt >overlap_raw_region.txt #fst_candidate_region.txt means region passed 1% threshold, B_ZM4_bed_region.txt means gene location, both file include chromosome, start and end columns
seqkit grep -f gene.list.txt  allgene.pep.fa >  candidate_gene.fa #use the candidate gene name extract protein info
##generate uniprot database, refer:https://www.uniprot.org/help/downloads
makeblastdb -dbtype prot -in  database/uniprot_sprot.fasta -out database/uniprot_sprot
blastp -query candidate_gene.fa -db database/uniprot_sprot -out candidate_gene.sprot.blast  -outfmt 6  -num_alignments 5  -num_threads 4
##please refer candidate_gene_filter.R for extract gene info
