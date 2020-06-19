# willow-workflows

#salix-ref-preparation
create index and collect preliminary statistics about reference

#willow_mapping
align sequence file to reference

#willow_preprocessing
collect preliminary statistics about samples bam file

#willow_call_mutation
call mutation sites using GATK

#snv_process
process to screen out candidate snv sites

#indel_process 
process to screen out candidate insertion and deletion sites

#screen-branch-specific-lack
pipeline to screen out snv sites(through the snv_process) shared by branches(BR-m)

#simulation
simulate mutation detect pipeline to evaluate this method's FNR
