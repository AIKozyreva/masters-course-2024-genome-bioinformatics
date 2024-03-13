# Repeats Masking in Genomes 

Сырые данные для работы: Melospiza melodia melodia isolate bMelMel2 chromosome Z - https://www.ncbi.nlm.nih.gov/nuccore/NC_086226.1/

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
Команда _conda install -c bioconda ncbi-tools-bin_ - don't work because of
```
Solving environment: failed
PackagesNotFoundError: The following packages are not available from current channels:
  - ncbi-tools-bin
Current channels:
  - https://conda.anaconda.org/bioconda
  - https://conda.anaconda.org/conda-forge
  - defaults
```
Итак, мы выяснили спустя 1.5 часа страданий, что инструкция просто старая, а пакета такого в conda не существует, и с похожим названием тоже не существует, 
найти его, может быть, получится как часть пакета NCBI c++ toolkit, с сайта: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/winmasker/README

Другой способ установить - это пойти в Google Colab или на локальную машину и установить через: (на сервере нужны права администратора, которых нет)
```
!sudo apt install ncbi-tools-bin
```
Но в колабе он работать тоже не будет, потому что после успешной установки вызов команды ncbi-tools-bin приведёт к следующему:
```
!sudo apt install ncbi-tools-bin
!install ncbi-tools-bin
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
_windowmasker [-h] [-help] [-xmlhelp] [-ustat unit_counts]
    [-in input_file_name] [-out output_file_name] [-checkdup check_duplicates]
    [-fa_list input_is_a_list] [-mem available_memory] [-unit unit_length] _
```

**<Рабочая схема установки на сервер dustmasker>**__
```
conda install blast=2.14.1
dustmasker -h                                                                                                                          
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
_Please use: trf File Match Mismatch Delta PM PI Minscore MaxPeriod [options]_
```
#### RepeatModeler and RepeatMasker
RepeatModeler and RepeatMasker can be installed together as they are often used in conjunction. Ensure to install dependencies like rmblast
```
conda install -c bioconda rmblast
conda install -c bioconda repeatmodeler repeatmasker
RepeatMasker -h
_RepeatMasker [-options] <seqfiles(s) in fasta format>_
RepeatModeler -h
_RepeatModeler [-options] -database <XDF Database>_
```
FINALLY!!!!!!!!!

### Step 3. Running Repeats Masking Tools

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
Output files are: trf_results.txt; NC_086226.1.fasta.2.7.7.80.10.50.500.dat   

#### RepeatModeler and RepeatMasker

Repeat Modeler - is a de novo transposable element (TE) family identification and modeling package. At the heart of RepeatModeler are three de-novo repeat finding programs ( RECON, RepeatScout and LtrHarvest/Ltr_retriever ) which employ complementary computational methods for identifying repeat element boundaries and family relationships from sequence data. 

Info: https://www.repeatmasker.org/RepeatModeler/ 
```
BuildDatabase -name vorobey_full_repeat_db -engine ncbi /mnt/projects/users/aalayeva/genomics/raw/GCF_035770615.1_bMelMel2.pri_genomic.fna
RepeatModeler -database vorobey_full_repeat_db -threads 16 -LTRStruct
_  <Using output directory = /mnt/projects/users/aalayeva/genomics/repeats/RM_2323907.TueMar121008052024
  Search Engine = rmblast 2.14.1+
  Threads = 16
  Dependencies: TRF 4.09, RECON , RepeatScout 1.0.6, RepeatMasker 4.1.5
  LTR Structural Analysis: Enabled (GenomeTools 1.6.4, LTR_Retriever, Ninja 0.97-cluster_only, MAFFT 7.520, CD-HIT 4.8.1 )
  Random Number Seed: 1710238085
  Database = /mnt/projects/users/aalayeva/genomics/repeats/vorobey_full_repeat_db .
  - Sequences = 353
  - Bases = 1541245192
  - N50 = 82773674>_
```
Repeat Masker - is a program that screens DNA sequences for interspersed repeats and low complexity DNA sequences. The output of the program is a detailed annotation of the repeats that are present in the query sequence as well as a modified version of the query sequence in which all the annotated repeats have been masked (default: replaced by Ns). There are two variants of command. Firsе one with custom database, made earlier by RepeatModeler, second one is with one of default databases, that was made for all birds.

Info: https://www.repeatmasker.org/webrepeatmaskerhelp.html 
```
1) RepeatMasker -lib /mnt/projects/users/aalayeva/genomics/repeats/RM_2323907.TueMar121008052024/consensi.fa -pa 4 -s -xsmall -e ncbi /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fa
2) RepeatMasker -species Birds -pa 4 -q -xsmall -e ncbi /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fa
```
Output files are: for variant 1) repeatmasker.table(1).txt, for var 2) repeatmasker.table.txt

For the second command i have got the better result: more places were maskered and programm identified more types of masked sequences. Probably, the reason of this were the termination of Repeat Modeler database building step. Or maybe the reason what in second variant Repeat Masker used default database, which includes information about many types of sequences of many Birds, while my custom database was based only on my one bird - Melospiza melodia - genome data.

### Step 4. Interpretation of Outputs

#### Table with output's fields comparision

| Tool/Fields | Window Masker | Tandem Repeats Finder | dust masker | Repeat Modeler | Repeat Masker |
| :---------- | :------------ | :-------------------- | :---------- | :------------- | :------------ |
|    Aim      | repeats masking |  repeats masking    | repeats masking | repeats identif-n | repeats masking |
|      1      |   True23.99   |     Codecademy Tee    |     False   |      23.99     |     23.99     |
|      2      |  False19.99   |     Codecademy Hoodie |     False   |      23.99     |     23.99     |
|      3      |  False42.99   |     Item              |    In Stock |      23.99     |     23.99     |




### Step 5. Discussion

Repeats masking is a crucial step in genome analysis pipelines due to the huge presence of repetitive elements within eukaryotic genomes. These repeats, such as transposons and retrotransposons or ShortTandemRep, can significantly compromise the accuracy of downstream analyses because failure to account for repetitive elements can lead to the algorithms false identification of multiple gene's copies or merge distinct genes due to shared repetitive sequences. This inaccuracy extends to regulatory regions, where unmasked repeats may corrupt the identification and interpretation of functional elements like promoters or enhancers. In evolutionary analyses, unmasked repeats can be misconstrued as conserved regions, affecting the delineation of homologous sequences and badly influent our understanding of species relationships.

