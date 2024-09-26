# Assembly
For the past year, I have been learning how to assemble a genome using ONT long-reads. At this moment, these scripts will serve to assemble a genome using different softwares, from correcting the reads, to scaffolding. This is an ON-GOING work, and therefore it has to be used with caution. I do not wish to make this an example pipeline, but more for organizing my work, and helping other peeps in the lab.

All the scripts here will be adequate for working with ONT long read data. I have no access to HIFi data as of now, so if you're working with that I recommend looking at other softwares. 

I have two possible pipelines in the scripts here:
1) Old pipeline:
- correct the reads using nanoq
- assemble using flye
- polish using medaka and pilon
- (opstional) phase using HapDup
- correct and scaffold using ragtag
- close gaps using TGS-Gap-Closer

2) New pipeline (will be named as "alternative")
- correct the reads using herro
- assemble using verkko (if you have ultra long reads) or hifiasm
- correct and scaffold using ragtag
- close gaps using TGS-Gap-Closer

I highly recommend using the new pipeline if you have the adequate data. This has given me best results for humans. The herro correction step is in the dorado new software and uses machine learning to correct the reads. It really works well, so I do not recommend doing any further polishing on the reads, especially if you're working with phased data. Polishing in later stages could really damage the phasing, since you're using non phased raw reads to polish your phased assembly. 

All of the links for specific githubs are in the scripts. At some point in the future I will add them all here. Please feel free to look at each script for now.
