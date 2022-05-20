# A transcriptome assembly from fragments of the annelids *Pygospio elegans* and *Arenicola marina* 

Authors:
* Elena Novikova
* Aleksandr Chen

## Introduction

*Annelids*, like many other invertebrate animals, replace lost body parts in a process called regeneration. In particular, *polychaetes*, a class of generally marine annelid worms, are capable of regenerating to some degree. The degree of regeneration varies widely across the taxon. For instance, *Pygospio elegans* (*Spionidae*, *Annelida*) is capable of regenerating both head and tail segments, whereas *Arenicola marina* (*Arenicolidae*, *Annelida*) does not regenerate lost segments. 

The concept behind *annelids’* ability to regenerate is called a positional information theory. This concept explains how cells determine their location in a multicellular structure. The positional information is conveyed by gradients of signaling molecules (proteins, RNA), that are produced in specific regions and induce concentration-dependent responses in target tissues.

## Aim and objectives

The **aim** is to assemble transcriptomes of two polychaetes - Pygospio elegans and Arenicola marina, - and to prepare data for further investigation of genes responsible for gradient expression in the body of annelids.

The following objectives were set in order to achive the goal:
1.  Data preparation 
    - Sequencing errors correction
    - Quality trimming and adapter clipping
    - Quality control 
2.  Transcriptome assembly
3.  Post-processing of assembly 
    - Reducing artificial sequence abundance
    - Decontamination
4.  Qualitative analysis of gene expression 
5.  Gene selection

## Data

Two samples (Pygospio elegans and Arenicola marina) were fragmented into 12 fragments. RNA from each of the fragment was sequnced. Overall, 24 sets of RNA-seq data or approximately 650 million of paired-end reads.

## Workflow
The workflow of the project presented below. Each part of scheme will be discussed later.
![Workflow](https://github.com/chensasha/TanscriptomeAssembly-Annelids/blob/main/Workflow.png)

---

## Data preparation
[Karect](https://github.com/aminallam/karect)
```
./karect -correct -threads=12 -matchtype=hamming -celltype=diploid -inputfile=PAIR_READS_1.fastq.gz -inputfile=PAIR_READS_2.fastq.gz
```

[Trimmomatic][http://www.usadellab.org/cms/?page=trimmomatic]
```
java -jar trimmomatic-0.39.jar PE CORR_PAIR_READS_1.fastq.gz CORR_PAIR_READS_2.fastq.gz output_forward_paired.fastq.gz output_forward_unpaired.fastq.gz output_reverse_paired.fastq.gz output_reverse_unpaired.fastq.gz ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:5:20 LEADING:25 TRAILING:25 MINLEN:25
```

## De novo transcriptome assembly

[Trinity][https://github.com/trinityrnaseq/trinityrnaseq/wiki]

```
Trinity --seqType fq --max_memory 100G --normalize_reads --left TRIMM_CORR_PAIR_READS_1.fastq.gz,... --right TRIMM_CORR_PAIR_READS_2.fastq.gz,... --CPU 6 --output PATH_TO_OUT_DIR
```

```
python rnaQUAST.py -c Trinity.fasta
```

## Post-processing of assembly 

[CD-HIT][http://weizhong-lab.ucsd.edu/cd-hit/]
```
cd-hit-est -i Trinity.fasta -o CDHIT -c 0.95 -M 2000
```

[barrnap][https://github.com/tseemann/barrnap]
```
barrnap --kingdom euk --threads 2 --outseq barrnap.fasta < CDHIT.fasta 
```

Get the uniref100 taxlist
```
git clone https://github.com/GDKO/uniref_taxlist.git
cat uniref100.taxlist.gz.part-a* | gunzip > uniref100.taxlist
```

Get and build uniref90 database with DIAMOND
```
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz | gunzip
diamond makedb --in uniref90.fasta --db uniref90
```

Download the ncbi taxonomy dmp
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | tar -xvf
```

```
bash MCSC_decontamination.sh file.ini
```

## Qualitative analysis of gene expression

```
salmon index -t CDHIT_clean.fasta -i CDHIT_clean_index -k 31 
```

```
salmon quant -i CDHIT_clean_index -l IU -1 TRIMM_CORR_PAIR_READS_1.fastq.gz ... -2 TRIMM_CORR_PAIR_READS_1.fastq.gz ... --validateMappings -o CDHIT_clean_salmon_out -p 8
```

```
./TransDecoder.LongOrfs -t CDHIT_clean.fasta
```

```
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
hmmsearch --cpu 8 --domtblout pfam.domtblout Pfam-A.hmm CDHIT_clean.fasta.transdecoder_dir/longest_orfs.pep
```

```
./TransDecoder.Predict -t CDHIT_clean.fasta --retain_pfam_hits pfam.domtblout 
```

## Gene selection
```
awk 'sub(/^>/, "")' Trinity.fasta > header.trinity # get the fasta headers
cat header.trinity | sed 's/_i.*//' > gene.header.trinity # keep only the portion of the header representing the gene_id
cat header.trinity | sed 's/[[:space:]]len.*//' > iso.header.trinity 
paste iso.header.trinity gene.header.trinity > Trinity.fasta.gene_trans_map.txt 
```

```
awk 'sub(/^>/, "")' CDHIT_Arenicola_clean.fasta.transdecoder.pep > headers_pep
cat headers_pep | sed 's/_i.*//' > gene_heades_pep
cat headers_pep | sed 's/.*len://' > pep_lengths
cat pep_lengths | sed 's/(.*//' > pep_lengths_all
paste gene_heades_pep pep_lengths_all > gene_pep_lengths
awk '{if ($2>=100) { print $1, $2 }}' gene_pep_lengths > pep_res
```
