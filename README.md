# Assembly
For the past year, I have been learning how to assemble a genome using ONT long-reads. At this moment, these scripts will serve to assemble a genome using different softwares, from correcting the reads, to scaffolding. This is an ON-GOING work, and therefore it has to be used with caution. 

All the scripts here will be adequate for working with ONT long read data. I have no access to HIFi data as of now, so if you're working with that I recommend looking at other softwares. 

I have two possible pipelines in the scripts here:
1) Old pipeline:
- correct the reads using nanoq
- assemble using flye
- polish using medaka and pilon
- correct and scaffold using ragtag

2) New pipeline (will be named as "alternative")
- correct the reads using herro
- assemble using verkko or hifiasm
- ON-GOING work

I highly recommend using the new pipeline if you have the adequate data. This has given me best results for humans, I have been able to resolve a very difficult genome region using hifiasm, and I am currently testing verkko.

All of the links for specific githubs are in the scripts. At some point in the future I will add them all here. Please feel free to look at each script for now.
