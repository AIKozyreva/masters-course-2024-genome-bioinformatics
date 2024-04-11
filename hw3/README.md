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
 господи прости я не могу, ну он не хочет
_Installation:_
```
////
```
_Run GeneMark on the genome:_
gms2.pl --seq <genome_file.fasta> --genome-type auto --output genemark_output.gff3
_Result:_

Predictions of genes are vital to validate against protein databases and transcriptomic alignments to ensure the accuracy and completeness of gene annotations. Comparing gene predictions with known protein sequences helps confirm the presence of coding regions and identify potential coding sequences missed during annotation. By aligning predicted genes with protein databases, researchers can assess whether predicted open reading frames correspond to known protein domains or motifs, providing confidence in gene structure and function assignments.

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
exonerate --query ./GCF_000211545.4_ASM21154v6_protein.faa --target annot.fna --model protein2genome --showalignment TRUE --ryo ">%ti %tab-%tae\n%tas\n" > Myc_exonerate_out_pr2gen.txt

```
_Result:_

> Query: WP_011113318.1 UMP kinase [Mycoplasmoides gallisepticum]
> Target: lcl|NC_023030.2 Mycoplasmoides gallisepticum chromosome, whole genome shotgun sequence
> Model: protein2genome:local
> Raw score: 1179
> Query range: 0 -> 236


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
mmseqs createdb ./annot.fna queryDB
mmseqs createdb ./GCF_000211545.4_ASM21154v6_protein.faa targetDB
```

_Running the search:_
```
mmseqs search proteinsDB genomeDB resultDB tmp
mmseqs convertalis proteinsDB genomeDB resultDB result.m8
head -n 10 resultDB.m8 
```
![image](https://github.com/AIKozyreva/masters-course-2024-genome-bioinformatics/assets/74992091/57b3020a-6ccd-4aab-9423-f77d59ba5cde)


By aligning annotated proteins with well-characterized ones, researchers can verify if annotations accurately represent functional domains, motifs, or features. This process helps in identifying potential errors or inaccuracies in the annotation, such as mislabeled domains or missed functional regions.
Moreover, protein alignment aids in detecting evolutionary relationships, highlighting conserved regions across species and indicating potential functional importance. Quality assessment through alignment also enables the identification of potential sequence errors, such as frameshifts or sequencing artifacts, which could affect downstream analyses. 

## Task 3: RNA-seq Mapping
Objective: Map RNA-seq data to the genome to provide transcript evidence for gene annotation.

### Alignment with Hisat2:
for mapping next-generation sequencing reads:
_Installation_
```
conda create -n hisat2
conda activate hisat2
conda install -c bioconda hisat2
```

_Build an index for your genome:_
```
hisat2-build <genome_file.fasta> genome_index
hisat2-build /mnt/projects/users/aalayeva/ref_human_38/GRch38.p14.genome.fa hisat2_GRCh38_INDEX
```
_Mapping RNA-seq reads to the genome:_
```
hisat2 -x genome_index -U <reads_file.fastq> -S mapped_reads.sam
hisat2 --fast --summary-file ./SUM.txt --reorder --qc-filter --add-chrname -k 1 -x ./index/hisat2_GRCh38_INDEX -U /mnt/projects/users/aalayeva/smallRNA/trim_10scamt.fastq -S ./align_hisat2.sam
```
Now we have to convert large .sam into binary .bam and sort it by cooordinates:
```
samtools view -bS align_hisat2.sam -o align_hisat2.bam
samtools sort align_hisat2.bam -o sorted_hisat2.bam
```

Transcriptomic alignments further validate gene predictions by aligning predicted transcripts with experimental RNA sequencing data, verifying their expression and splicing patterns. This comparison ensures that predicted genes accurately represent transcriptional activity under different conditions and developmental stages. Moreover, transcriptomic alignments help detect alternative splicing events or untranslated regions missed during gene prediction, enhancing the completeness of gene annotations.

Overall, cross-referencing gene predictions with protein databases and transcriptomic alignments is essential for refining gene annotations, improving their accuracy, and enhancing our understanding of gene structure, function, and regulation.

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
