#!/bin/sh
#   brad-ref-preparation.sh - Preparing of reference dataset.
#
#   Author: Nowind
#   Created: 2017-04-13
#   Updated: 2017-04-13
#
#   References:
#   http://www.broadinstitute.org/gatk
#   http://bio-bwa.sourceforge.net/
#   http://samtools.sourceforge.net/
#
#   Change logs:
#   Version 1.0.0 17/04/13: The initial version.
#



echo ===============[$0 v$VERSION]==============
echo `date` ": start process..."
echo



######################################################################
## build indexes

## bwa index
bwa index /ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa 2>&1 | \
    tee /ref/Spurpurea/assembly/Spurpurea_289_v1.0.bwa-index.log


## samtools faidx
samtools faidx /ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa


## picard CreateSequenceDictionary
java -jar /opt/nfs/share/biosoft/picard-tools-1.114/CreateSequenceDictionary.jar \
    REFERENCE=/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
    OUTPUT=/ref/Spurpurea/assembly/Spurpurea_289_v1.0.dict

######################################################################




######################################################################
## nucleotide content


## nucleotide number, GC content
parseSEQs.pl --fasta /ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
    > /ref/Spurpurea/assembly/Spurpurea_289_v1.0.content.csv

parseSEQs.pl --sum-all --fasta /ref/Spurpurea/assembly/Spurpurea_289_v1.0.softmasked.fa \
    > /ref/Spurpurea/assembly/Spurpurea_289_v1.0.softmasked.content.csv


## triplets
count_triplets.pl --fasta /ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
    > /ref/Spurpurea/assembly/Spurpurea_289_v1.0.triplets.csv
######################################################################




