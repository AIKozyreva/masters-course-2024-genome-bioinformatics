Homework 3. Gene Annotation in Genomes
This assignment aims to provide you with practical experience in gene annotation processes, combining ab initio prediction, homology-based approaches, and transcript evidence. By completing this assignment, you will learn how to:

Predict genes in a genome using AUGUSTUS and GeneMark.
Use Exonerate and MMSeqs2 for protein alignment to find homologous proteins.
Map RNA-seq data to the genome for transcript evidence.
Compare the annotated genes with RefSeq proteins to validate and refine annotations.

References:
Genes annotation part 1 colab: https://colab.research.google.com/drive/1CZM2J8gHKpJuXcNjZBFqvTvUh4yu63V4
Genes annotation part 2 colab: https://colab.research.google.com/drive/1a2sWhDVDwME1YGojHEoDqDT2tUw-oEgM
Assignment: https://github.com/aglabx/masters-course-2024-genome-bioinformatics/blob/main/homeworks/hw3_2024.md 

## Task 1: Gene Prediction 
Objective: Use AUGUSTUS and GeneMark for ab initio gene prediction in the provided genome sequence.

### AUGUSTUS:
_Installation:_
```
conda create -n augustus
conda activate augustus
conda install -c bioconda augustus
```

_Run AUGUSTUS on the provided genome:_
```
augustus --species=chicken --gff3=on --progress=true /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fasta > ./genomics/augustus_output.gff3
```
_Result:_
```
(base) aalayeva@aglab1:~/genomics$ head -n 10 augustus_output.gff3 
> ##gff-version 3
> #This output was generated with AUGUSTUS (version 3.3.3).
```

### GeneMark:
эту хрень невозможно установить из-за ключей, кто вообще написал всё на перле, а потом ещё и намутил с лицензиями господи прости 

_Installation:_
```

```

_Run GeneMark on the genome:_

gms2.pl --seq <genome_file.fasta> --genome-type auto --output genemark_output.gff3

_Result:_


## Task 2: Homology-Based Annotation 
Objective: Use Exonerate and MMSeqs2 to align proteins from a closely related species to the genome.

### EXONERATE 

_Installation_
```
conda create -n exonerate
conda activate exonerate
conda install -c bioconda exonerate
```

_Align proteins to the genome:_
```
exonerate --model protein2genome --showtargetgff True <protein_sequences.fasta> <genome_file.fasta> > exonerate_output.gff3
```
_Result:_



### MMSeqs2 for Proteins:

![ ..](https://github.com/AIKozyreva/masters-course-2024-genome-bioinformatics/blob/main/hw3/%D0%91%D0%B5%D0%B7%D1%8B%D0%BC%D1%8F%D0%BD%D0%BD%D1%8B%D0%B9.png) 

_Installation:_
```
conda create -n mmseqs2
conda activate mmseqs2
conda install -c bioconda mmseqs2
```

_Creating a database for the protein sequences and the genome:_
```
mmseqs createdb <protein_sequences.fasta> proteinsDB
mmseqs createdb <genome_file.fasta> genomeDB
mmseqs search proteinsDB genomeDB resultDB tmp
```

_Running the search:_
```
mmseqs convertalis proteinsDB genomeDB resultDB result.m8
```

Include a discussion on the significance of the protein alignments in gene annotation.

## Task 3: RNA-seq Mapping
Objective: Map RNA-seq data to the genome to provide transcript evidence for gene annotation.

### Alignment with Hisat2:
for mapping next-generation sequencing reads:
_Installation_
```

conda install -c bioconda hisat2
```

_Build an index for your genome:_
```
hisat2-build <genome_file.fasta> genome_index

```
_Mapping RNA-seq reads to the genome:_
```
hisat2 -x genome_index -U <reads_file.fastq> -S mapped_reads.sam
hisat2 --fast --summary-file ./SUM.txt --reorder --qc-filter --add-chrname -k 1 -x ./index/hisat2_GRCh38_INDEX -U /mnt/projects/users/aalayeva/smallRNA/trim_10scamt.fastq -S ./align_hisat2.sam
```

Discuss how RNA-seq evidence supports gene annotation.

## Task 4: Comparing Annotations with RefSeq GFF.
**Tools: BEDtools, python, seqkit, bioawk**

Objective: Validate and refine your gene annotations by comparing the intersection of intervals between your annotation results and RefSeq GFF data.
Instructions:
This task requires you to compare your annotated genes (in GFF format) with the RefSeq GFF file for your organism of interest. 
The goal is to understand how well your annotations match up with the RefSeq annotations, which serves as a standard reference.

_Bedtools Installation:_

Install Bedtools, a powerful tool for genome arithmetic, which you will use to compare GFF files:
```
conda install -c bioconda bedtools
```

Intersecting GFF Files:
Use Bedtools to find overlaps between your GFF file and the RefSeq GFF file.
This step will help identify matching annotations, partial overlaps, and unique annotations in both files.

bedtools intersect -a predicted_genes.gff3 -b refseq_genes.gff3 -wa -wb > intersection_results.txt

Analysis:
Open the resulting intersection_results.txt file.
Analyze the intersections to evaluate the coverage and accuracy of your annotations compared to RefSeq.
Identify and document cases of exact match, partial overlap, and unique annotations in your report.

Deliverables:
The report should include the command used for intersecting GFF files.
A summary of the intersection analysis:
The number of exact matches, partial overlaps, and unique annotations.
Discussion on the significance of these results for the quality of your annotation.
Additional Notes:
For detailed instructions and examples on how to perform these tasks, see Genes annotation part 1 colab. This Colab notebook provides step-by-step guidance on handling GFF files, executing intersection operations, and analyzing the results effectively. Make sure to adapt the instructions from the Colab to fit the specific details of your assignment, such as file paths and organism names.

Note to remember:

As you work through this task, consider the implications of the intersections for your annotations' accuracy and completeness. Document any discrepancies between your annotations and the RefSeq data, as these can highlight areas for further investigation or refinement in your gene prediction methods.

Good luck with your comparative analysis, and may it enrich your understanding of gene annotation accuracy and reliability!

Expected files in hw3 folder:
README.md
Predicted gene files (*.gff3)
Alignment and mapping results .bam

Report file (README.md) should include:
Conda commands used for environment setup and tool installations.
Commands used for gene prediction, protein alignment, RNA-seq mapping, and comparison with RefSeq proteins.
Summary and analysis of each task's results.
