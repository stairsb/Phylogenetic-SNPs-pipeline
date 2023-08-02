# Phylogenetic-SNPs-pipeline
This pipeline pulls out SNPs from whole genome data and results in a phylogenetic tree
## Aligning reads to Regerence Genome
Aligned reads to reference genome using: https://bio-bwa.sourceforge.net/bwa.shtml
```
bwa meme /path/to/reference/genome/.fasta /path/to/read/one.fastq /path/to/read/two/fastq > output.sam
```
## Creating human readible and binary files
Convert the `sam` to `bam` which will convert the file to binary format to allow for faster computing. Index and sort the `bam` which will make for access specific short reads faster. http://www.htslib.org/doc/samtools.html
```
samtools view -bS output.sam > output.bam
samtools sort output.bam > output_sorted.bam
samtools index output_sorted.bam
```
## SNP Calling
We will use `bcftools mpileup` to generate a file conatining the genotype likelihoods. `bcftools call` is what we will use to call the SNPs. Finally, `bcftools view` can be used to view, subset and filter VCF or BCF files by position and filtering expression. 
Download `vsfutils.pl` found here https://github.com/lh3/samtools/blob/master/bcftools/vcfutils.pl. bcftools documentation: http://www.htslib.org/doc/bcftools.html#view
```
bcftools mpileup -f pathtoreferencegenome.fasta pathto_sorted.bam -o isolate_paired_sorted.vcf
bcftools call --ploidy 1 -c pathto_paired_sorted.vcf -o isolate_paired_sorted.bcf
bcftools view isolate_paired_sorted.bcf | pathto/vsfutils.pl varFilter -D 200 > outfile_SNPoutput
```
## Merging SNPs
The following script can be found in the GWAS scripts repository. Merge all of the SNPs into a single file:
```
MergeSNPs.py SNPoutput.fasta > SNP.txt
```
Filter the alleles into different files based on their allele frequencies. 
```
Biallelic_SNPs.py -o merged -m SNP.txt
```
To get all the SNPs into fasta format run. Input files name should be `merged_SNPs.txt`
```
SNP_fasta.py merged_SNPs.txt
```
## Phylogenetic Analysis
We will be using iqtree to perform our phylogenetic analysis. iqtree documentation: http://www.iqtree.org/
```
iqtree -b 1000 -m GTR+ASC+F -nt 8 -s SNPs.fasta
```
## Plotting in R-studio
R packges needed
```
library(ggtree)
library(phytools)
library(ggnewscale)
library(ggrepel)
```

Read in phylogenetic reconstruction and metadata
```
SNP_tree2 <- read.tree("micro_SNPs.fasta.treefile")
data2 <- read.csv("micro_RM_names.csv")
```
Format tip labels and assign metadata
```
SNP_tree2$tip.label<-gsub("_"," ",SNP_tree2$tip.label)
group_names2 <- list(Clinical = c(), Environmental = c(), Unknown = c())

for (row in 1:nrow(data2)) {
  if(data2[row, "Isolation.Site"] == "Clinical"){
    print(data2[row, "Isolation.Site"])
    group_names2[["Clinical"]] <- append(group_names2[["Clinical"]], data2[row, "Strain"]) 
  }  
  else if (data2[row, "Isolation.Site"] == "Environmental"){
    group_names2[["Environmental"]] <- append(group_names2[["Environmental"]], data2[row, "Strain"])
  }
  else {
     group_names2[["Unknown"]] <- append(group_names2[["Unknown"]], data2[row, "Strain"])
  }
}
```

Reroot tree to the right outgroup
```
rerooted <- reroot(SNP_tree2, 23)
rerooted$node.label <- gsub("Root"," ",rerooted$node.label)
grouping_ecology <- groupOTU(rerooted, group_names2, group_name = "Ecology")
```
Add some color to make things look how we want
```
pal <- wes_palette("Cavalcanti1", 5, type = "continuous")
```

Actual phylogenetic tree
```
SNP_p <- ggtree(grouping_ecology) +
  geom_hilight(node = 87, fill = "#000000", alpha = 0.1, extend = 0.025) +
  geom_highlight(node = 88, fill = "#333333", alpha = 0.1, extend = 0.025) +
  geom_highlight(node = 89, fill = "#666666", alpha = 0.05, extend = 0.21) +
  scale_color_manual(values=c("red", "green4", "black"), labels = c("Clinical", "Environmental"), na.translate=FALSE) +
  theme(legend.key.size = 10) +
  geom_tree() +
  theme_tree() +
  geom_treescale() +
  geom_rootedge(0.03) +
  guides(color = guide_legend(override.aes = list(label = "\u25CF", size = 5))) +
  geom_tiplab(size = 3, aes(color=Ecology)) +
  xlim(0,0.28) +
  new_scale_color()

 SNP_p2 <- SNP_p %<+% data2 +
  geom_tippoint(aes(colour=Variety, size = Shape)) + 
  scale_size_continuous(range = c(1.5)) +
  scale_colour_discrete(type = pal, limits = c("azygosporus", "chinensis", "microsporus", "oligosporus", "rhizopodiformis")) +
  guides(size = FALSE) +
  theme(legend.position = c(0.896, 0.836), legend.spacing.y = unit(0, "cm"))

SNP_p2
```
Finally, the output file is saved. There should be multiple output file formats
```
ggsave("SNP.pdf", device = cairo_pdf, width = 30, height = 20, units = "cm", limitsize = FALSE)
```



