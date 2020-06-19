**willow-workflows**
==================================================
**prepare and mapping**
---------------------------------------------
### salix-ref-preparation
shell script to create index and collect preliminary statistics about reference

### willow_mapping
shell script to align sequence file to reference

### willow_preprocessing
shell script to collect preliminary statistics about samples bam file

**detecr mutation**
--------------------------------------------
### willow_call_mutation
shell script to call mutation sites using GATK

### snv_process
shell script to process to screen out candidate snv sites

### indel_process 
shell script to process to screen out candidate insertion and deletion sites

### screen-branch-specific-lack
perl script to screen out snv sites(through the snv_process) shared by branches(BR-m)

**simulation**
---------------------------------------------------
### simulation
shell script to simulate mutation detect pipeline to evaluate this method's FNR
