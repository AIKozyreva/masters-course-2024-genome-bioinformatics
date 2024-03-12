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
```
sed -i 's/\[.*\]//g' /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fa
windowmasker -in /mnt/projects/users/aalayeva/genomics/raw/GCA_035770615.1.fa -mk_counts -out ./windowMasker_genome.counts
windowmasker -in /mnt/projects/users/aalayeva/genomics/raw/GCA_035770615.1.fa -ustat ./windowMasker_genome.counts -out ./windowmasker_results.txt -outfmt fasta
```
Output files are: windowmasker_results.txt and windowMasker_genome.counts

```
dustmasker -in /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fasta -out /mnt/projects/users/aalayeva/genomics/repeats/dustmasker_result.fasta -outfmt fasta -parse_seqids
dustmasker -in /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fasta -out /mnt/projects/users/aalayeva/genomics/repeats/dustmasker_res_infoasn1.txt -outfmt maskinfo_asn1_text
```
Output files are: dustmasker_res_infoasn1.txt and dustmasker_result.fasta

#### Tandem Repeats Finder (TRF)
```
trf /mnt/projects/users/aalayeva/genomics/raw/NC_086226.1.fasta 2 7 7 80 10 50 500 -d -h > /mnt/projects/users/aalayeva/genomics/repeats/trf_results.txt
```
Output files are: trf_results.txt; NC_086226.1.fasta.2.7.7.80.10.50.500.dat   

#### RepeatModeler and RepeatMasker
```


```


