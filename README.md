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

