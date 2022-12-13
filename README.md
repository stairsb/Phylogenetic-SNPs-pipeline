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
bcftools view pathto/vsfutils.pl varFilter -D 200 > outfile_SNPoutput
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
Read in a tree in newick format and define the groups for each isolate that is included in the tree.
```
stree <- read.tree("SNPs.fasta.treefile")
cls2 <- list(Clinical=c("B07675","B07643","B07367","B06600","B06590","B05459","B07675","B07643","B07585","B07386","B10187","B08956","B11541","B11543","B11546","B11547","B11535","B11538","B11539","B11540","B11553","B11550","B11554","B10881","B10548","B11523","B11147","B11529","B11526","B11532","B11531","B11534","B11533","B11552","B11551","B11557","B11556","B11555"), Environmental=c("5550","5547","13129","5546","5547","5548","5550","5551","5552","5553","5554","5558","17693","B11549"))
stree2 <- groupOTU(stree, cls2)
```
Tree with branches all visiable. The branch distances between taxon are not based on caluated values but this tree does help to show the relationship between taxon.
```
ggtree(stree2, branch.length = "none") + 
  geom_tiplab(size = 1.8, aes(color=group))  +  
  theme(legend.position="right") + 
  geom_nodepoint(size = 0.6) + 
  scale_color_manual(values=c("red", "green4")) + 
  geom_strip('B11555', 'B11538', barsize = 1, color='chocolate',     label = "Argentina", offset.text = .1, fontsize = 2, offset = 0.6) + 
  geom_strip('17693', '5554', barsize = 1, color='blue3', label = "USA", offset.text = .1, fontsize = 2, offset = 0.6) + 
  geom_strip('B11533', 'B11534', barsize = 1, color='hotpink2', label = "Argentina", offset.text   = .1, fontsize = 2, offset = 0.6) +
  geom_strip('B10548', 'B07643', barsize = 1, color='blue3', label = "USA", offset.text   = .1, fontsize = 2, offset = 0.6) +
  geom_strip('B11147', 'B11543', barsize = 1, color='hotpink2', label = "Argentina", offset.text   = .1, fontsize = 2, offset = 0.6) +
  geom_strip('B10187', 'B07367', barsize = 1, color='blue3', label = "USA", offset.text   = .1, fontsize = 2, offset = 0.6) +
  geom_strip('B11541', 'B11557', barsize = 1, color='hotpink2', label = "Argentina", offset.text   = .1, fontsize = 2, offset = 0.6) +
  geom_strip('5550', 'B07585', barsize = 1, color='blue3', label = "USA", offset.text   = .1, fontsize = 2, offset = 0.6) +
  geom_strip('5547', '5547', barsize = 1, color='darkgoldenrod', label = "Brazil", offset.text = .1, fontsize = 2, offset = 0.6) +      
  geom_strip('5546', '5546', barsize = 1, color='chocolate4', label = "Philippenes", offset.text =    .1, fontsize = 2, offset =  0.6) +
  geom_strip('B07386', '13129', barsize = 1, color='blue3', label = "USA", offset.text   = .1, fontsize = 2, offset = 0.6) + 
  geom_text2(aes(subset = !isTip, label=label), size = 1.8, nudge_x = .28)
```
![image](https://user-images.githubusercontent.com/111078377/207306066-d2b51718-1632-42dd-93d1-5e997bfc3484.png)

Actual tree
```
ggtree(stree2) + 
  geom_tiplab(align = TRUE, linesize = 0.5, size = 1.8, aes(color=group))  + 
  theme(legend.position="right") + 
  geom_nodepoint(size = 0.6) + 
  scale_color_manual(values=c("red", "green4")) + 
  geom_strip('B11555', 'B11538', barsize = 1, color='hotpink2', label = "Argentina", offset.text = .01, fontsize = 2, offset = 0.009) + 
  geom_strip('17693', '5554', barsize = 1, color='blue3', label = "USA", offset.text = .01, fontsize = 2, offset = 0.009) + 
  geom_strip('B11533', 'B11534', barsize = 1, color='hotpink2', label = "Argentina", offset.text   = .01, fontsize = 2, offset = 0.009) +
  geom_strip('B10548', 'B07643', barsize = 1, color='blue3', label = "USA", offset.text   = .01, fontsize = 2, offset = 0.009) +
  geom_strip('B11147', 'B11543', barsize = 1, color='hotpink2', label = "Argentina", offset.text   = .01, fontsize = 2, offset = 0.009) +
  geom_strip('B10187', 'B07367', barsize = 1, color='blue3', label = "USA", offset.text   = .01, fontsize = 2, offset = 0.009) +
  geom_strip('B11541', 'B11557', barsize = 1, color='hotpink2', label = "Argentina", offset.text   = .01, fontsize = 2, offset = 0.009) +
  geom_strip('5550', 'B07585', barsize = 1, color='blue3', label = "USA", offset.text   = .01, fontsize = 2, offset = 0.009) +
  geom_strip('5547', '5547', barsize = 1, color='darkgoldenrod', label = "Brazil", offset.text = .01, fontsize = 2, offset = 0.009) +      
  geom_strip('5546', '5546', barsize = 1, color='chocolate4', label = "Philippenes", offset.text =    .01, fontsize = 2, offset = 0.009) +
  geom_strip('B07386', '13129', barsize = 1, color='blue3', label = "USA", offset.text   = .01, fontsize = 2, offset = 0.009)
```
![image](https://user-images.githubusercontent.com/111078377/207306252-dcec9103-f77e-4030-9a9c-56a1296964c2.png)




