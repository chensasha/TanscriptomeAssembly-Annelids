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
Error correction is done using [Karect](https://github.com/aminallam/karect).
```
./karect -correct -threads=12 -matchtype=hamming -celltype=diploid -inputfile=PAIR_READS_1.fastq.gz -inputfile=PAIR_READS_2.fastq.gz
```

Then, low-quality and adapter sequences are clipped
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
```
java -jar trimmomatic-0.39.jar PE CORR_PAIR_READS_1.fastq.gz CORR_PAIR_READS_2.fastq.gz output_forward_paired.fastq.gz output_forward_unpaired.fastq.gz output_reverse_paired.fastq.gz output_reverse_unpaired.fastq.gz ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:5:20 LEADING:25 TRAILING:25 MINLEN:25
```

## De novo transcriptome assembly

The reads are now ready for de novo assembly. They are assembled via [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) using all prepared libraries.
```
Trinity --seqType fq --max_memory 100G --normalize_reads --left TRIMM_CORR_PAIR_READS_1.fastq.gz,... --right TRIMM_CORR_PAIR_READS_2.fastq.gz,... --CPU 6 --output PATH_TO_OUT_DIR
```

We can now check some statistics of the assembly with [rnaQUAST](https://github.com/ablab/rnaquast) to get a quality assessment.
```
python rnaQUAST.py -c Trinity.fasta
```

Some statistics are shown below.
![Nx](https://github.com/chensasha/TanscriptomeAssembly-Annelids/blob/main/rnaQUAST/Nx.png)
![Length](https://github.com/chensasha/TanscriptomeAssembly-Annelids/blob/main/rnaQUAST/transcript_length.png)

## Post-processing of assembly 

Trinity de novo assembly has artificial abundance due to use of De Bruijn graphs. In order to get rid of that we use [CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/) that clusters similar sequneces into clusters.
```
cd-hit-est -i Trinity.fasta -o CDHIT -c 0.95 -M 2000
```

Then, we predict the location of ribosomal RNA genes using [barrnap](https://github.com/tseemann/barrnap).
```
barrnap --kingdom euk --threads 2 --outseq barrnap.fasta < CDHIT.fasta 
```
NCBI search revealed contaminants like *Selenidium pygospionis*. So decontamination needs to be conducted. We're going to use [MCSC](https://github.com/Lafond-LapalmeJ/MCSC_Decontamination). 

Before we can run MCSC, we need to get the uniref100 taxlist.
```
git clone https://github.com/GDKO/uniref_taxlist.git
cat uniref100.taxlist.gz.part-a* | gunzip > uniref100.taxlist
```

Then, get and build uniref90 database with [DIAMOND](https://github.com/bbuchfink/diamond).
```
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz | gunzip
diamond makedb --in uniref90.fasta --db uniref90
```

Download the ncbi taxonomy dmp.
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | tar -xvf
```

Now we can run MCSC with the parameters listed in file.ini (setting taxon to keep as 'Annelida').
```
bash MCSC_decontamination.sh file.ini
```

## Qualitative analysis of gene expression

Now that we have clustered decontaminated data we can analyze gene expression. We use [Salmon](https://github.com/COMBINE-lab/salmon) to produce transcript-level quantification estimates from our data. First, build an index.
```
salmon index -t CDHIT_clean.fasta -i CDHIT_clean_index -k 31 
```

Then, we can quantify our set our reads.
```
salmon quant -i CDHIT_clean_index -l IU -1 TRIMM_CORR_PAIR_READS_1.fastq.gz ... -2 TRIMM_CORR_PAIR_READS_1.fastq.gz ... --validateMappings -o CDHIT_clean_salmon_out -p 8
```

We now identify candidate coding regions within transcript sequence via [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki). First, extract the long open reading frames.
```
./TransDecoder.LongOrfs -t CDHIT_clean.fasta
```

Search the peptides for protein domains using [Pfam](http://pfam.xfam.org) and [hmmer3](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan).
```
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
hmmsearch --cpu 8 --domtblout pfam.domtblout Pfam-A.hmm CDHIT_clean.fasta.transdecoder_dir/longest_orfs.pep
```

We can now predict the likely coding regions. The outputs generated above can be leveraged by TransDecoder to ensure that those peptides with domain hits are retained in the set of reported likely coding regions. 
```
./TransDecoder.Predict -t CDHIT_clean.fasta --retain_pfam_hits pfam.domtblout 
```

## Gene selection
We can finally select genes with significant expression. 

First, we need to get a file of Trinity transcripts and corresponding genes.
```
awk 'sub(/^>/, "")' Trinity.fasta > header.trinity # get the fasta headers
cat header.trinity | sed 's/_i.*//' > gene.header.trinity # keep only the portion of the header representing the gene_id
cat header.trinity | sed 's/[[:space:]]len.*//' > iso.header.trinity # keep the transcript name
paste iso.header.trinity gene.header.trinity > Trinity.fasta.gene_trans_map.txt # combine two columns
```
Now we use a library [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) for R to summarize expression to genes by rinning Expression_sum.R script. Out of those genes we select those who have >1 TPM expression.

Next, we choose genes that encode proteins longer than 100 amino acids.
```
awk 'sub(/^>/, "")' CDHIT_clean.fasta.transdecoder.pep > headers_pep # get the fasta headers
cat headers_pep | sed 's/_i.*//' > gene_heades_pep # keep only the portion of the header representing the gene_id
cat headers_pep | sed 's/.*len://' > pep_lengths # keep the lengths of the proteins
cat pep_lengths | sed 's/(.*//' > pep_lengths_all 
paste gene_heades_pep pep_lengths_all > gene_pep_lengths # combine two columns
awk '{if ($2>=100) { print $1, $2 }}' gene_pep_lengths > pep_res 
```

## Results
We assembled two transcriptomes from 12 fragments of *Pygospio elegans* and *Arenicola marina*. Boreover, two sets of protein-coding genes selected with significant expression: 54315 (Pygospio) and 33530 (Arenicola).

## Future plans
* Prepared data, Pygospio elagans and Arenicola marina assemblies and sets of genes, can be further analized in order to determine gene-candidates responsible for positional information concept.

* Moroever, there are two more samples of considered annelids available with sequnecing data. The same process can be apllied in order to validate the results.

## Literature
1.  Heikkinen LK, Kesäniemi JE, Knott KE. De novo transcriptome assembly and developmental mode specific gene expression of Pygospio elegans. Evol Dev. 2017 Jul;19(4-5):205-217.
2.  Joël Lafond-Lapalme,  Marc-Olivier Duceppe,  Shengrui Wang,  Peter Moffett, Benjamin Mimee. A new method for decontamination of de novo transcriptomes using a hierarchical clustering algorithm. Bioinformatics, Volume 33, Issue 9, 1 May 2017, Pages 1293–1300.
