**willow-workflows**
==================================================
**prepare and mapping**
---------------------------------------------
### 1.salix-ref-preparation.sh
shell script to create index and collect preliminary statistics about reference

### 2.willow_mapping.sh
shell script to align sequence file to reference

### 3.willow_preprocessing.sh
shell script to collect preliminary statistics about samples bam file

**detect mutation**
--------------------------------------------
### 1.willow_call_mutation.sh
shell script to call mutation sites using GATK

### 2.snv_process.sh
shell script to screen out candidate snv sites

### 3.indel_process.sh 
shell script to process to screen out candidate insertion and deletion sites

### 4.screen-branch-specific-lack.pl
perl script to screen out snv sites(through the snv_process) shared by branches(BR-m)

**simulation**
---------------------------------------------------
### 1.simulation.sh
shell script to simulate mutation detect pipeline to evaluate this method's FNR
