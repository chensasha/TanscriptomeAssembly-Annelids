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

![Workflow](https://github.com/chensasha/TanscriptomeAssembly-Annelids/blob/main/Workflow.png)

