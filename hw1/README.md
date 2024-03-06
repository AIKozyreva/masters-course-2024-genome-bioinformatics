#### Requirements for selected genome:
- refseq annotation NCBI number - GCF_035770615.1 
- Selected genome - Melospiza melodia melodia
- Певчая овсянка — вид певчих воробьиных птиц из семейства Passerellidae, живущих в Северной Америке. one of the most abundant, variable and adaptable species.
- NCBI Genome: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_035770615.1/
- USCS database: https://genome.ucsc.edu/cgi-bin/hgGateway?hgsid=1980412654_uIe6IvhMHaMO0lMtpXArvHzAabjP&redirect=manual&source=genome.ucsc.edu
- Ensembl db: absence

#### General info about in_files: out_gen_info.txt
Amount of files lines: 353 and 354 lines, headers has a different info, ncbi files headers contain information about samples and sequencing type, uscs headers looks like id's (maybe chromosomes id's).

#### Quast report - report.html
```
conda install -c bioconda quast
quast.py --eukaryote --large --fast --silent /content/GCA_035770615.1.fa
```
##### QUAST Results:
      
|            Metrics                  |   Value    | 
| ----------------------------------- | ---------- |
|            Largest contig bp        | 174870981  |
|            Total length             | 1541261971 |
|                 N50                 | 77460261   |
|                 L50                 |     7      | 
|                 GC (%)              |   44.96    |

  - Largest contig is the longest whole sequence
  - N50 это длина последовательности, совокупность последовательностей такой или большей длины точно покрывает 50% (total length) всей последовательности собранных контигов.
  - L50 is the number of contigs equal to or longer than N50. L50 - is the minimal number of contigs that cover half the assembly.

#### BUSCO report 
```
conda create -n busco -c conda-forge -c bioconda busco=5.6.1
tar -xzf /mnt/projects/users/aalayeva/genomics/busco/busco_downloads/passeriformes_odb10.2024-01-08.tar.gz
busco -i /mnt/projects/users/aalayeva/genomics/raw/GCA_035770615.1.fa -m genome -l /mnt/projects/users/aalayeva/genomics/busco/passeriformes_odb10 --offline --download_path /mnt/projects/users/aalayeva/genomics/busco/busco_downloads
```
##### Results:
      
|            Metrics                  |   Value  | 
| ----------------------------------- | -------- |
|                Complete BUSCOs (C)  | 10401    |
| Complete and single-copy BUSCOs (S) | 10343    |
| Complete and duplicated BUSCOs (D)  |   58     |
|            Fragmented BUSCOs (F)    |   74     | 
|               Missing BUSCOs (M)    |   369    |
|      Total BUSCO groups searched    |   10844  | 

  - Complete BUSCOs (C) Genes that are present, complete, and not fragmented
  - Fragmented BUSCOs (F) Genes that are present but not complete
  - Missing BUSCOs (M) Genes are expected to be present based on the reference, but not

##### BUSCO Assembly Statistics:

|            Metrics                  |   Value    | 
| ----------------------------------- | ---------- |
|            Number of scaffolds      |    354     |
|            Total length bp          | 1541261971 |
|            Percent gaps             |   0.070%   |
|            Scaffold N50             |   77 MB    | 

