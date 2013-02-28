## Introduction

In this exercise we try to find all SNPs. These single-nucleotide-polymorphisms are the main cause for variantion in genomes. SNPs in genes are interesting, because they are the reason for variants or mutations of a gene.

## 7.2.1 Detect the SNPs (SNP calling) using SAMtools and the same reference genome

1. Use SAMtools to generate a VCF (Variant Call Format) file.

		The samtools man-page provides us with the following procedure to collect the data:	

		samtools mpileup -ugf Trinity.fasta Illumina.sorted.q20.bam | bcftools view -bvcg - > var.raw.bcf
		bcftools view var.raw.bcf | vcfutils.pl varFilter -d5 -D1000 > var.raw.vcf

2. Explore the .vcf file. What columns would you be interested in?
		
		For this exercised the very first column is sufficient; it gives me the names of the genes, which have a SNP.
		
3. Use grep to count the number of SNPs. How many SNPs have been called?
		
		The command `grep -c '^comp' var.raw.vcf` shows, that 162286 SNPs have been found.

4. Count the number of SNPs per transcript. Create a table with two columns (contig-ID and count) and name it SNP_counts.txt. Hint: A combination of grep, cut, sort, and uniq could be useful.  
		
		The command `egrep -o "comp[0-9]+\_c0\_seq[1-2]" var.raw.vcf | uniq -c | sed 's/^[ ]*//g' | awk '{print $2"\t"$1}' > SNP_counts.txt` creates the file with the gene-ids and SNP-counts.

5. Import the table you generated (i.e., the number of SNPs per transcript) into R and merge with your expression data.

		>snp_counts = read.delim("SNP_counts.txt", header=FALSE, sep="")
		>merged_data <- read.csv2("merged_data.csv")
		>merged_data = merge (merged_data, snp_counts, by = "V1")
		>colnames(merged_data)[15] = "SNP.counts"
		>merged_data$SNPdensity = (merged_data$SNP.count /
		>merged_data$length)*1000


6. Explore your data.

		Commands in R and their output:		

		>names(merged_data)
		 [1] "V1"         "X"          "counts.454" "counts.Ill" "Annotation"
		 [6] "attribute"  "start"      "end"        "score"      "strand"    
		[11] "frame"      "gene_ID"    "length"     "FPKM"       "SNP.counts"
		[16] "SNPdensity"

		>head (merged_data)
		                 V1  X counts.454 counts.Ill Annotation attribute start  end
		1 comp10012_c0_seq1  1          0       4896    Ensembl      exon   100 1844
		2 comp10030_c0_seq1  2          0        235    Ensembl      exon   100 1922
		3 comp10035_c0_seq1  3          0        244    Ensembl      exon   100 2782
		4 comp10041_c0_seq1  4          0        491    Ensembl      exon   100 2601
		5 comp10070_c0_seq1 10          0        306    Ensembl      exon   100 1858
		6 comp10084_c0_seq1 12          0       4155    Ensembl      exon   100 2450
		  score strand frame                   gene_ID length       FPKM SNP.counts
		1    42      +     . gene_id=comp10012_c0_seq1   1744 22837.7863          2
		2    42      +     . gene_id=comp10030_c0_seq1   1822  1049.2490          2
		3    42      +     . gene_id=comp10035_c0_seq1   2682   740.0995          8
		4    42      +     . gene_id=comp10041_c0_seq1   2501  1597.0808          6
		5    42      +     . gene_id=comp10070_c0_seq1   1758  1415.9947          2
		6    42      +     . gene_id=comp10084_c0_seq1   2350 14383.4227         21
		  SNPdensity
		1   1.146789
		2   1.097695
		3   2.982849
		4   2.399040
		5   1.137656
		6   8.936170


7. How was the SNP density calculated?
		The R-code "merged_data$SNPdensity = (merged_data$SNP.count / merged_data$length)*1000" shows, that the formula is: (SNP.counts * 1000) / length

8. Get the gene with highest density of SNPs.

		>head (merged_data[order(merged_data$SNPdensity, decreasing
		+ = TRUE),])
		                   V1   X counts.454 counts.Ill Annotation attribute start  end
		19  comp10232_c0_seq1  37          0        274    Ensembl      exon   100 2253
		133 comp11975_c0_seq2 545          0       2934    Ensembl      exon   100 1768
		44  comp10670_c0_seq1 114          0        254    Ensembl      exon   100 1882
		112 comp11777_c0_seq1 416          0      32594    Ensembl      exon   100 2138
		103 comp11642_c0_seq1 355          0        370    Ensembl      exon   100 2704
		113 comp11777_c0_seq2 417          0      17780    Ensembl      exon   100 2123
		    score strand frame                   gene_ID length       FPKM SNP.counts
		19     42      +     . gene_id=comp10232_c0_seq1   2153   1035.299         60
		133    42      +     . gene_id=comp11975_c0_seq2   1668  14309.457         38
		44     42      +     . gene_id=comp10670_c0_seq1   1782   1159.538         36
		112    42      +     . gene_id=comp11777_c0_seq1   2038 130104.572         34
		103    42      +     . gene_id=comp11642_c0_seq1   2604   1155.899         41
		113    42      +     . gene_id=comp11777_c0_seq2   2023  71498.174         31
		    SNPdensity
		19    27.86809
		133   22.78177
		44    20.20202
		112   16.68302
		103   15.74501
		113   15.32378

		The gene comp10232_c0_seq1 has the highest SNPdensity.

9. Humans have an average of 1 SNP per kb while Drosophila has about 1 SNP every 100 bp. What is the average nucleotide diversity (#SNPs/kb) in the O. vulgaris transcriptome?

		>summary (merged_data$SNPdensity)
		   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
		 0.1865  0.8106  1.6550  2.6790  3.1020 27.8700 

10. Would you say it is more similar to humans or to Drosophila?

		We calculate the SNPdensity for humans and drosophila:
		* human: 1 SNP per 1kbp: 			1 * 1000 / 1000 	= 1
		* drosophila: 1 SNP per 100bp: 	1 * 1000 / 100 	= 10

		The median and the mean from the summary are lying much more close to the snp.count to humans than to that of drosophila.

## Discussion
In terms of SMPs our genome of the octopus vulgaris is more closely related to mammals than to insects.