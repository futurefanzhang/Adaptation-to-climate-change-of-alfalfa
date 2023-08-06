
mydata=read.table("candidate_gene.sprot.blast",header = F,sep = "\t")
colnames(mydata)=c("qseqid", "sseqid",    "pident",    "length",    "mismatch",    "gapopen",    "qstart",    "qend",    "sstart",    "send",    "evalue",    "bitscore")
#filter p-value<1e-10
mydata=mydata[mydata$evalue<1e-10,]
mydata=mydata[order(mydata$pident,decreasing = TRUE),]
mydata=mydata[!duplicated(mydata$qseqid),] #remove duplicated gene info
write.table(mydata,"deg.sprot.blast_use.txt",sep = "\t",quote = F,row.names = F,col.names = F) #export use blast info
write.table(mydata$sseqid,"gene_function_and_Uniprot_info.txt",sep = "\t",quote = F,row.names = F) #export UniProt accessions info, use DAVID web check function,https://david.ncifcrf.gov/list.jsp
