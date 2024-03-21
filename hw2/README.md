# Repeats Masking in Genomes 

Сырые данные для работы: Melospiza melodia melodia isolate bMelMel2 chromosome Z - https://www.ncbi.nlm.nih.gov/nuccore/NC_086226.1/

## Task 0: Theory about genome-repeats

**Tandem Repeats:**
Definition: Tandem repeats are sequences of DNA where identical or nearly identical nucleotide motifs are repeated successively and appear adjacent to each other along the DNA strand.
Organization: The repeated units are typically clustered together in a head-to-tail fashion.
Examples: Microsatellites or simple sequence repeats (SSRs), minisatellites, and variable number tandem repeats (VNTRs) are common examples of tandem repeats.

**Interspersed Repeats:**
Definition: Interspersed repeats are sequences that are scattered throughout the genome rather than being clustered together.
Organization: Unlike tandem repeats, interspersed repeats are not organized in a head-to-tail fashion but are dispersed at various locations.
Examples: Transposons, retrotransposons, and other mobile genetic elements fall under interspersed repeats.

**Low-Complexity Sequences:**
Definition: Low-complexity sequences are regions in DNA that consist of simple and repetitive motifs but may not necessarily be organized in tandem or interspersed patterns.
Organization: These repeats can occur in various arrangements, and their organization may not follow a specific pattern.
Examples: Regions with homopolymeric runs, such as poly-A or poly-T stretches, are examples of low-complexity sequences.

**Satellite DNA:**
Definition: Satellite DNA consists of highly repetitive sequences that are often organized into longer arrays, forming distinct bands during centrifugation.
Organization: Satellite DNA can be tandemly arranged (forming satellites) or interspersed in certain regions.
Examples: Centromeric and telomeric repeats are types of satellite DNA.

## Task 1: Installing Repeats Masking Tools. 
Objective: Install WindowMasker, DUST, TRF, RepeatModeler, and RepeatMasker using Conda to familiarize yourself with software installation and environment management.

### Step1. Conda Environment Setup
```
conda create -n repeatmasker
conda activate repeatmasker
conda install python=3.8
```
### Step 2. Installing Tools
#### WindowMasker and DUST
WindowMasker and DUST are part of the NCBI toolkit. 
Команда _conda install -c bioconda ncbi-tools-bin_ - **don't work because** of
```
conda install -c bioconda ncbi-tools-bin
____________________________________________________________________
Solving environment: failed
PackagesNotFoundError: The following packages are not available from current channels:
  - ncbi-tools-bin
```
Итак, мы выяснили спустя 1.5 часа страданий, что инструкция просто старая, а пакета такого в conda не существует, и с похожим названием тоже не существует, 
найти его, может быть, получится как часть пакета NCBI c++ toolkit, с сайта: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/winmasker/README

Другой способ установить - это пойти в Google Colab или на локальную машину: (на сервере нужны права администратора, которых нет)
```
!sudo apt install ncbi-tools-bin
```
Но в колабе **он работать тоже не будет**, потому что после успешной установки вызов команды ncbi-tools-bin приведёт к следующему:
```
!sudo apt install ncbi-tools-bin
!install ncbi-tools-bin
_______________________________________________________________________
_/bin/bash: line 1: ncbi-tools-bin: command not found_
```
**<Так в итоге рабочая схема установки на сервер WindowMasker>**__
1. Создаём директорию, куда положим windowmasker и переходим в неё
2. Скачиваем прямой ссылкой
3. Делаем файл исполняемым
4. Вызывеем программу
```
mkdir ./windowmasker
cd ./windowmasker
wget ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/windowmasker/windowmasker 
chmod +x windowmasker
./windowmasker -h
______________________________________________________________________________
_windowmasker [-h] [-help] [-xmlhelp] [-ustat unit_counts]
    [-in input_file_name] [-out output_file_name] [-checkdup check_duplicates]
    [-fa_list input_is_a_list] [-mem available_memory] [-unit unit_length] _
```

**<Рабочая схема установки на сервер dustmasker>**__
```
conda install blast=2.14.1
dustmasker -h
_______________________________________________________________________
_dustmasker [-h] [-help] [-xmlhelp] [-in input_file_name]                                                                                                               
    [-out output_file_name] [-window window_size] [-level level]                                                                                                 
    [-linker linker] [-infmt input_format] [-outfmt output_format]                                                                                                       
    [-parse_seqids] [-hard_masking] [-version-full] [-version-full-xml]                                                                                                  
    [-version-full-json] _                                                                                                                                                                                 
```

#### Tandem Repeats Finder (TRF)
```
conda install -c bioconda trf
trf -h
___________________________________________
_Please use: trf File Match Mismatch Delta PM PI Minscore MaxPeriod [options]_
```
#### RepeatModeler and RepeatMasker
RepeatModeler and RepeatMasker can be installed together as they are often used in conjunction. Ensure to install dependencies like rmblast
```
conda install -c bioconda rmblast
conda install -c bioconda repeatmodeler repeatmasker
RepeatMasker -h
RepeatModeler -h
_______________________________________________________
_RepeatMasker [-options] <seqfiles(s) in fasta format>_
_______________________________________________________
_RepeatModeler [-options] -database <XDF Database>_
```
_**FINALLY!!!!!!!!!**_

### Task 3. Running Repeats Masking Tools

#### WindowMasker and DUST
WindowMasker - is a program that identifies and masks out highly repetitive DNA sequences and DNA sequences with low complexity in a genome using only the sequence of the genome itself.

Info: https://github.com/hydrachallenge/main/blob/master/manuals/README.windowmasker 
```
sed -i 's/\[.*\]//g' /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fa
windowmasker -in /mnt/projects/users/aalayeva/genomics/raw/GCA_035770615.1.fa -mk_counts -out ./windowMasker_genome.counts
windowmasker -in /mnt/projects/users/aalayeva/genomics/raw/GCA_035770615.1.fa -ustat ./windowMasker_genome.counts -out ./windowmasker_results.txt -outfmt fasta
```
Output files are: windowmasker_results.txt and windowMasker_genome.counts

DustMasker is a program that identifies and masks out low complexity parts of a genome using a new and improved DUST algorithm. The main advantages of the new algorithm are symmetry with respect to taking reverse complements, context insensitivity, and much better performance

Info: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/README 
```
dustmasker -in /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fasta -out /mnt/projects/users/aalayeva/genomics/repeats/dustmasker_result.fasta -outfmt fasta -parse_seqids
dustmasker -in /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fasta -out /mnt/projects/users/aalayeva/genomics/repeats/dustmasker_res_infoasn1.txt -outfmt maskinfo_asn1_text
```
Output files are: dustmasker_res_infoasn1.txt and dustmasker_result.fasta

#### Tandem Repeats Finder (TRF)
TRF - is a program to locate and display tandem repeats in DNA sequences. In order to use the program, the user submits a sequence in FASTA format. There is no need to specify the pattern, the size of the pattern or any other parameter. The output consists of two files: a repeat table file and an alignment file. The repeat table contains information about each repeat, including its location, size, number of copies and nucleotide content. Repeats with pattern size in the range from 1 to 2000 bases are detected.

Info: https://tandem.bu.edu/trf/definitions
```
trf /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fasta 2 7 7 80 10 50 500 -d -h > /mnt/projects/users/aalayeva/genomics/repeats/trf_results.txt
```
Output files are: trf_results.txt; NC_086226.1.fasta.2.7.7.80.10.50.500.dat ; and also can be a lot of .html (if there was no -h parametr ib command). 
Summary .html can looks like this one and it's the same information as .dat file consists of, but in sum.html each link connects with other smaller .html, where you will find information about detected repeats in one or another sequnce (contig, chromosome or something other)

![TFR-html-output](https://github.com/AIKozyreva/masters-course-2024-genome-bioinformatics/blob/main/hw2/TRF_html_out_example.jpg?raw=true)

#### RepeatModeler and RepeatMasker

Repeat Modeler - is a de novo transposable element (TE) family identification and modeling package. At the heart of RepeatModeler are three de-novo repeat finding programs ( RECON, RepeatScout and LtrHarvest/Ltr_retriever ) which employ complementary computational methods for identifying repeat element boundaries and family relationships from sequence data. 

Info: https://www.repeatmasker.org/RepeatModeler/ 
```
BuildDatabase -name vorobey_full_repeat_db -engine ncbi /mnt/projects/users/aalayeva/genomics/raw/GCF_035770615.1_bMelMel2.pri_genomic.fna
RepeatModeler -database vorobey_full_repeat_db -threads 16 -LTRStruct
_____________________________________________________________________________
_<Search Engine = rmblast 2.14.1+
Dependencies: TRF 4.09, RECON , RepeatScout 1.0.6, RepeatMasker 4.1.5
LTR Structural Analysis: Enabled (GenomeTools 1.6.4, LTR_Retriever, Ninja 0.97-cluster_only, MAFFT 7.520, CD-HIT 4.8.1 )
Database = /mnt/projects/users/aalayeva/genomics/repeats/vorobey_full_repeat_db>_
```
Repeat Masker - is a program that screens DNA sequences for interspersed repeats and low complexity DNA sequences. The output of the program is a detailed annotation of the repeats that are present in the query sequence as well as a modified version of the query sequence in which all the annotated repeats have been masked (default: replaced by Ns). There are two variants of command. Firsе one with custom database, made earlier by RepeatModeler, second one is with one of default databases, that was made for all birds.

Info: https://www.repeatmasker.org/webrepeatmaskerhelp.html 
```
1) RepeatMasker -lib /mnt/projects/users/aalayeva/genomics/repeats/RM_2323907.TueMar121008052024/consensi.fa -pa 4 -s -xsmall -e ncbi /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fa

2) RepeatMasker -species Birds -pa 4 -q -xsmall -e ncbi /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fa
```
Output files are: for variant 1) repeatmasker.table(1).txt, for var 2) repeatmasker.table.txt

For the second command i have got the better result: despite on more places were maskered, programm identified more types of masked sequences -as i see that's more reliable than first variant, where programm masked more sequences about 10% (based on lib, got from Modeler round4), but all masked sequences were Unclassified or Simple repeats, while with command with default db we can see a lot of differen types of maskered sequences.

Probably, the reason of this were the termination of Repeat Modeler second step after round4. Or maybe the reason - is that in second variant Repeat Masker has used default database, which includes information about many types of sequences of many Birds, while my custom database was based only on my one bird - Melospiza melodia - genome data.

### Task 4. Interpretation of Outputs

#### dustmasker
- _'fasta' - genome sequences with masked subseq-s_
- _'interval' - genes and all masked groups in them. only by coordinates_
- _'maskinfo_asn1_text' - list of masked groups with their iDs, start end coordinates in bp_
```
>lcl|NC_086226.1_gene_1 [gene=LOC134431993] [db_xref=GeneID:134431993] [location=complement(14..2329)] [gbkey=Gene]
266 - 272
509 - 516
591 - 599
```
You can create a sctipt which will count the masked groups for each gene and their summary, but there is no information about groups type. 

```
awk 'BEGIN { gene_count = 0; total_count = 0; } /^>/ { if (gene_count > 0) { print gene_count; gene_count = 0; } print; next } { gene_count++; total_count++; } END { print gene_count; print "General amount of masked groups is", total_count, "for", NR, "genes"; }' /mnt/projects/users/aalayeva/genomics/repeats/dustmasker_result_1303.txt > /mnt/projects/users/aalayeva/genomics/repeats/out_sum_dustmasker

tail -n 3  /mnt/projects/users/aalayeva/genomics/repeats/out_sum_dustmasker
```
**General amount of masked groups is 53260 for 54232 sequences**

#### windowmasker
- '.counts' is the output of WindowMasker Stage 1 processing. It also serves as input for Stage 2 processing. The first line of the file contains one integer number which is the unit size. Then come the lines containing counts for the units which appeared more than T_low times in the genome (and its reverse complement). These are ordered by the unit numerical value.
- '-outfmt interval' the output of stage2 WindowMasker is Interval format. The file is consisting of blocks of information for each input sequence in the input FASTA order.  Each block starts with the FASTA title of the seq followed by the description of masked intervals, one interval per line. The intervals do not overlap and are sorted by their start position. (NOTE: the positions are numbered starting at 0.)

- _'windowMasker_genome.counts' - list of counts for the units which appeared more than T_low times in the genome_
```
##pct: 500 100
12
0 68
```
- _'windowmasker_results.txt' - genome sequences with masked subseq-s_
```
>lcl|NC_086226.1_gene_3
acttgttttaaaattctacTGTAAGCGTTTATggttcactttttttttctttgcagcacc
cgtgacttttttttccctctaaaTCGGGAGTAAATATAGCTGAGTAAAATTAccttgatt
ttcttccttcttttttaactgaattgtttttattaaaaatcctaTTTGAATATGTAGAAA
```

#### Tandem Repeats Finder (TRF)

_'.dat' - The summary table includes the following information the pattern of repeat and the whole masked subsequence in this point (coordinates in the begining) of this gene_

```
Sequence: lcl|NC_086226.1_gene_8 
Parameters: 2 7 7 80 10 50 500
1406 1451 25 1.9 25 86 8 69 58 4 4 32 1.37 ATAAATTGATACAATATAAAATAAC ATAAATTGTAAATATAAAATAACATAAATTGATACAATATATAATA
2081 2133 6 9.3 6 73 19 51 52 0 47 0 1.00 AAAGGG AAAGGGGAAGGGAAAGGAAAGGGAAGGAAAGGGAAAAGGAAAAGGGAAGGGAA
```
#### RepeatModeler

- _'families.fa' - Consensus sequences_
- _'<database_name>-families.stk' - Seed alignments_
- _'<database_name>-rmod.log' - summarized log of the run_

```
RepeatModeler Round # 5
========================
Searching for Repeats
 -- Sampling from the database...
   - Gathering up to 270000000 bp
   - Sequence extraction : 00:10:20 (hh:mm:ss) Elapsed Time
 -- Running TRFMask on the sequence...
   - TRFMask time 00:28:36 (hh:mm:ss) Elapsed Time
 -- Masking repeats from the previous rounds...
       118916 repeats masked totaling 59597810 bp(s).
   - TE Masking time 00:14:16 (hh:mm:ss) Elapsed Time
 -- Sample Stats:
       Sample Size 270310507 bp
       Num Contigs Represented = 187
       Non ambiguous bp:
             Initial: 270002066 bp
             After Masking: 171227757 bp
             Masked: 36.58 %
 -- Input Database Coverage: 440479992 bp out of 1541245192 bp ( 28.58 % )
Sampling Time: 00:53:30 (hh:mm:ss) Elapsed Time
Running all-by-other comparisons...
  - Total Comparisons = 23028291
```

#### RepeatMasker

- _'RepeatMasker.fa.masked'-ordinary fasta masked sequences_

- _'RepeatMasker.fa.tbl'-main summary file_
```
file name: NC_086226.1.fa
sequences:           972
total length:   37852161 bp  (37839288 bp excl N/X-runs)
GC level:         39.90 %
bases masked:    2933633 bp ( 7.75 %)
```
![RepeatMasker.fa.tbl - output_example](https://github.com/AIKozyreva/masters-course-2024-genome-bioinformatics/blob/main/hw2/repmask_out_tbl_example.jpg?raw=true)

![RepeatMasker.out - output_example](https://github.com/AIKozyreva/masters-course-2024-genome-bioinformatics/blob/main/hw2/repeatmasker_out_example.jpg?raw=true)

- _'RepeatMasker.ori.out' - file with coord, gene and pattern of repeat_
```
23 26.01 1.54 0.00 lcl|NC_086226.1_gene_1      1633      1697    (619) +               (ATA)n   Simple_repeat       1      66     (0)
15 13.16 0.00 5.71 lcl|NC_086226.1_gene_10       396       432   (6024) +              (TTTG)n   Simple_repeat       1      35     (0)
```
#### Table with output's fields comparision
 
| Tool/Fields |  Window Masker  | Tandem Repeats Finder |   dust masker   |     Repeat Modeler    |  Repeat Masker  |
| :---------- | :-------------- | :-------------------- | :-------------- | :-------------------- | :-------------- |
|    Aim      | repeats masking |    repeats masking    | repeats masking |   repeats identif-n   | repeats id-n + masking |
| input_frmt  |      Fasta      |         Fasta         |      Fasta      |    Fasta + database   | Fasta      |
| output_frmt |        txt      |       .dat = text     |  Fasta & other  | .fasta + .str + .log  | fa.masked+ori.out+.tbl+.out |
| columns     | coord of masked groups/gene | pattern, sequence, coord of repeat | coord of masked groups/gene | db + repeats_seq-s | pattern, type, sequence, coord, stats of repeat |

| Tool/Fields |  Window Masker  | Tandem Repeats Finder |   dust masker   |  Repeat Masker  |
| :---------- | :-------------- | :------------------- | :-------------- | :-------------- |
| Total bp    |    37852161     |       37852161       |    37852161     |    37852161     |
| Masked bp   |     4926663     |         521458       |     1299929     |      3866474    |
| Masked %    |     13.02%      |         1.38%        |     3.43%       |      10.21%     |
| Mask  type  |    lowercase    |      N-masking       |     lowercase   |    lowercase    |
| gap bp      |      12873      |         none         |       none      |       12873     |
| GAP %       |      0.03%      |         none         |       none      |       0.03%     |
### Task 5. Discussion

**Why we have to mask something?**
Repeats masking is a crucial step in genome analysis pipelines due to the huge presence of repetitive elements within eukaryotic genomes. These repeats, such as transposons and retrotransposons or ShortTandemRep, can significantly compromise the accuracy of downstream analyses because failure to account for repetitive elements can lead to the algorithms false identification of multiple gene's copies or merge distinct genes due to shared repetitive sequences. This inaccuracy extends to regulatory regions, where unmasked repeats may corrupt the identification and interpretation of functional elements like promoters or enhancers. In evolutionary analyses, unmasked repeats can be misconstrued as conserved regions, affecting the delineation of homologous sequences and badly influent our understanding of species relationships.

- **WindowMasker**
Pros:
Efficient at identifying low-complexity regions and masking them.
Suitable for quickly identifying and masking repetitive elements.
Cons:
May lack the sensitivity to detect more complex repetitive structures.
Limited functionality compared to comprehensive repeat analysis tools.

- **DustMasker**
Pros:
Efficient at identifying and masking low-complexity regions.
Quick and computationally less intensive.
Cons:
Primarily designed for identifying low-complexity sequences, may not capture more diverse repetitive elements.
Limited in its ability to identify and mask interspersed repeats.

- **Tandem Repeat Finder**
Pros:
Specialized in identifying tandem repeats, making it effective for certain genomic regions.
Provides information on the pattern and organization of tandem repeats.
Cons:
Focuses solely on tandem repeats and may not capture interspersed repeats or complex repetitive structures.
Sensitivity can be affected by parameter settings. HARDMASKING - replacing repat sequences by NNNN make the process of detecting gaps in assembly more harder. Also hardmasking leads to the information losing.

- **RepeatModeler**
Pros:
De novo identification and modeling of repetitive elements in a genome.
Integrates multiple algorithms for a more comprehensive analysis.
Suitable for identifying novel repeats.
Cons:
Computationally more intensive compared to simple masking tools.
Requires more computational resources for de novo repeat identification.

- **RepeatMasker**
Pros:
Integrates RepeatModeler's repeat libraries for masking.
Provides detailed reports on identified repeats and their locations.
Customizable to include user-defined repeat libraries.
Cons:
Can be resource-intensive for large genomes.


In terms of algorithmic relationships, RepeatModeler often uses de novo algorithms for repeat identification, and RepeatMasker utilizes the generated repeat libraries for masking. Windowmasker and dustmasker are the easiest tools, while others arre more comlex and includes more algorithms, which works together (based on windowmasker or dust alg-s too). //maybe, actually i'm not sure, as i understood//
