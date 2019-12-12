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
bwa index /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa 2>&1 | \
    tee /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.bwa-index.log


## samtools faidx
samtools faidx /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa


## picard CreateSequenceDictionary
java -jar /opt/nfs/share/biosoft/picard-tools-1.114/CreateSequenceDictionary.jar \
    REFERENCE=/home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
    OUTPUT=/home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.dict


## blastdb
makeblastdb -dbtype nucl -in /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
    -title Bdistachyon_314_v3.0 \
    -out /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa
######################################################################




######################################################################
## nucleotide content


## nucleotide number, GC content
parseSEQs.pl --fasta /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
    > /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.content.csv

parseSEQs.pl --sum-all --fasta /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.softmasked.fa \
    > /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.softmasked.content.csv


## triplets
count_triplets.pl --fasta /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
    > /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.triplets.csv
######################################################################




######################################################################
## repeat content

## masked infos
awk 'BEGIN{OFS="\t";} !/\#/ {print $2,$4,$5;}' /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0_centromere.csv \
    > /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0_centromere.bed

maskedSEQ2bed.pl /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
    > /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.base_N.bed

cat /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0_centromere.bed \
    /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.base_N.bed | \
    sort -k1,1 -k2,3n | bedtools merge -i \
    > /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0_centromere_with_N.bed



## Tandem Repeats Finder
## http://tandem.bu.edu/trf/trf.unix.help.html
## The recomended values for Match Mismatch and Delta are 2, 7, and 7 respectively;
## The best performance can be achieved with values of PM=80 and PI=10
time trf /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa 2 7 7 80 10 20 150 -m -f -d -h

convert_trf.pl -i /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa.2.7.7.80.10.20.150.dat \
    -f bed | sort -k1,1 -k2,3n | bedtools merge -i - -nms \
    > /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.trf.w150.bed
######################################################################



######################################################################
## gene annotations

gff2tables.pl --gff /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.gff3 \
    --all-features --no-funcs | sed 's/_v2\.0\.a1//g' \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.tbl

awk 'BEGIN{OFS="\t";} !/\#/ && $1 ~ /Pp/ && $2 ~ /CDS/ {print $1,($3-1),$4,$5;}' \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.tbl | sort -k1,1 -k2,3n | \
    bedtools merge -i - -nms \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.cds.bed

awk '{print $1"\t"$3-$2; print "Total\t"$3-$2;}' /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.cds.bed | \
    subtotal_stats.pl -i - \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.cds.stats.csv

perl -ne 'next if (/\#/ || /scaffold/); my ($chrom, $source, $type, $start, $end, $strand, $info) = (split /\s+/)[0..4,6,8];
    next unless($type eq "CDS"); $info =~ /ID\=(.*?)(.v2.1.CDS.\d+\;|$)/; my $name = $1; $start -= 1;
    print "$chrom\t$start\t$end\t$name\n";' /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.gene.gff3 | \
    sort -k1,1 -k2,3n | bedtools merge -i - -nms \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.gene.cds.bed

awk '{print $1"\t"$3-$2; print "Total\t"$3-$2;}' /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.gene.cds.bed | \
    subtotal_stats.pl -i - \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.gene.cds.stats.csv


parseSEQs.pl --fasta /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds_primaryTranscriptOnly.fa \
    --count-codon --skip "Prupe.I" \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds_primaryTranscriptOnly.codon_counts.csv

parseSEQs.pl --fasta /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds_primaryTranscriptOnly.fa \
    --rate --skip "Prupe.I" --sum-all \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds_primaryTranscriptOnly.contents.csv


## extract promoter regions
awk 'BEGIN{OFS="\t";} !/\#/ && $1 ~ /Pp/ && $2 ~ /promoter/ {print $1,($3-1),$4,$2,$7;}' \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.tbl \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.promoter.bed


## extract UTR regions
awk 'BEGIN{OFS="\t";} !/\#/ && $1 ~ /Pp/ && $2 ~ /UTR/ {print $1,($3-1),$4,$2,$7;}' \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.tbl \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.UTR.bed


## extract gene regions (without UTRs)
awk 'BEGIN{OFS="\t";} !/\#/ && $1 ~ /Pp/ && $2 ~ /mRNA/ {print $1,($3-1),$4,$2,$7;}' \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.tbl \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.mRNA.bed

cat /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.promoter.bed \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.UTR.bed | \
    bedtools subtract -a /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.mRNA.bed -b - | sed 's/mRNA/TR/g' \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.TR.bed


## extract intergenic regions
cat /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.mRNA.bed \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.promoter.bed \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.UTR.bed | \
    bedtools complement -i - -g /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.chrom.fai | \
    awk 'BEGIN{OFS="\t";} {print $1,$2,$3,"intergenic\tinter"NR;}' \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.intergenic.bed

cat /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.promoter.bed \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.UTR.bed \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.TR.bed \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.intergenic.bed \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.features.bed

awk '{print $4"\t"($3-$2);}' /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.features.bed | \
    subtotal_stats.pl -i - > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.features.stats.csv


cat /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.promoter.bed \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.UTR.bed \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.TR.bed | \
    awk 'BEGIN{OFS="\t";} {$2 += 1; print $1,$2,$3,$5":"$4;}' | sort -k1,1 -k2,3n | bedtools merge -n -nms -i - \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/genes/Prunus_persica_v2.0.a1.primaryTrs.features.olp.bed




##
## format function annotations
##
perl -e 'my %funcs = (); print "#locusName\tDescriptions\n";
    while(<>){next if (/^\#/); chomp; my ($locusName, $Id, $IdType, $description) = (split /\t/);
    $description ||= "-"; push @{$funcs{$locusName}->{$IdType}}, "$Id,$description";}
    for my $locusName (sort keys %funcs) {
        my @annotations = ();
        for my $IdType (qw(PFAM PANTHER GO KEGGORTH KOG EC)) {
            my $descs = $funcs{$locusName}->{$IdType} ? (join "|", sort @{$funcs{$locusName}->{$IdType}}) : "-";
               $descs =~ s/\s+/_/g; push @annotations, "$IdType=$descs";
        } my $annotations = join ";", @annotations; print "$locusName\t$annotations\n";}' \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/functional/Prunus_persica_v2.0.a1_gene_functions.txt \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/functional/Prunus_persica_v2.0.a1_gene_functions.csv


perl -ne 'next if (/\#/); my ($chrom, $source, $type, $start, $end, $strand, $info) = (split /\s+/)[0..4,6,8];
    next unless($type eq "gene"); $info =~ /Name\=(.*?)(\;|$)/; my $length = $end - $start + 1;
    my $name = $1; print "$source\t$chrom\t$type\t$start\t$end\t$length\t$name\t$strand\n";' \
    /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.gene.gff3 | \
    map_records.pl --subject /home/wl/Data/data/peach/ref/GDR_v2.0.a1/functional/Prunus_persica_v2.0.a1_gene_functions.csv \
    --query - --rows1 6 --rows2 0 | \
    awk 'BEGIN{OFS="\t"; print "#source\tchrom\ttype\tstart\tend\tlength\tgene\tstrand\tdescriptions";}
    {if($9 == "N/A"){$10 = "NO_Predicted_Functions";}print $1,$2,$3,$4,$5,$6,$7,$8,$10;}' \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.gene.csv

awk '!/\#/' /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.gene.csv | \
    cut -f 1,3 --complement \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.gene.bed


## generate custom annotation files for BINGO
## http://www.psb.ugent.be/cbd/papers/BiNGO/Customize.html
perl -e 'print "(species=Prunus persica)(type=Full)(curator=GDR)\n";
    while(<>){next if (/^\#/); chomp; my ($locusName, $Id, $IdType, $description) = (split /\t/);
    next unless($IdType eq "GO"); $Id =~ /GO:(\d+)/; my $go_number = $1; print "$locusName = $go_number\n";}' \
    /home/wl/Data/data/peach/ref/GDR_v2.0.a1/functional/Prunus_persica_v2.0.a1_gene_functions.txt | \
    awk 'NR == 1; NR > 1 {print $0 | "sort"}' \
    > /home/wl/Data/data/peach/ref/GDR_v2.0.a1/functional/Prunus_persica_v2.0.a1_gene_functions.BINGO.dat





##
## local pfam
##
hmmscan --noali --cpu 2 \
    --tblout /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.tbl \
    --domtblout /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.domtbl \
    --pfamtblout /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.pfamtbl \
    -o /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam \
    /home/wl/Data/data/datasets/Pfam/Pfam28.0/Pfam-A.hmm \
    /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.fa

perl -e 'my %funcs = (); print "#locusName\tPFAM\n";
    while(<>){next if (/^\#/); chomp; my @line = (split /\s+/);
    my $description = join "\_", @line[18..$#line]; push @{$funcs{$line[2]}}, "$line[0]:$description";}
    for my $locusName (sort keys %funcs) {
        my $pfam_out = join ";", (sort @{$funcs{$locusName}}); print "$locusName\t$pfam_out\n";}' \
    /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.tbl \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.csv



##
## local interpro
##
sed 's/\*$//' /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.fa \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.fasta

time interproscan --input /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.fasta \
    --formats TSV --disable-precalc --seqtype p --tempdir /home/wl/Data/tmp/ \
    --outfile /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.interpro.tsv

perl -e 'my %funcs = (); print "#locusName\tDescriptions\n";
    while(<>){next if (/^\#/); chomp; my ($locusName, $IdType, $Id, $description) = (split /\t/)[0,3..5];
    $description ||= "-"; push @{$funcs{$locusName}->{$IdType}}, "$Id,$description";}
    for my $locusName (sort keys %funcs) {
        my @annotations = ();
        for my $IdType (qw(Pfam ProSiteProfiles ProSitePatterns Hamap SMART TIGRFAM)) {
            my $descs = $funcs{$locusName}->{$IdType} ? (join ";", sort @{$funcs{$locusName}->{$IdType}}) : "-";
               $descs =~ s/\s+/_/g; push @annotations, "$IdType=$descs";
        } my $annotations = join ";", @annotations; print "$locusName\t$annotations\n";}' \
    /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.interpro.tsv \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.interpro.csv



##
## F-box
##
awk '/F-box/ {print $3;}' \
    /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.tbl | sort | uniq \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.F-box.txt

fasta_process.pl --query /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.F-box.txt \
    --fasta /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds.fa --rows 0 | \
    parseSEQs.pl --fasta - \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds.fa.F-box.length.csv



##
## LRR-Pkinase, LRR-only
##
awk '/LRR/ {print $3;}' \
    /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.tbl | sort | uniq \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.LRR.txt

awk '!/#/ && !/LRR/ {print $3;}' \
    /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.tbl | sort | uniq \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.non-LRR.txt

map_records.pl --query /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.LRR.txt \
    --subject /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.non-LRR.txt \
    --rows1 0 --rows2 0 | awk '$2 == "N/A" {print $1;}' \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.LRR-only.txt

fasta_process.pl --query /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.LRR-only.txt \
    --fasta /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds.fa --rows 0 | \
    parseSEQs.pl --fasta - \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds.LRR-only.length.csv


awk '/Pkinase/ {print $3;}' \
    /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.tbl | sort | uniq \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.Pkinase.txt

map_records.pl --query /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.LRR.txt \
    --subject /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.Pkinase.txt \
    --rows1 0 --rows2 0 | awk '$2 != "N/A" {print $1;}' \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.LRR-Pkinase.txt


fasta_process.pl --query /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.LRR-Pkinase.txt \
    --fasta /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds.fa --rows 0 | \
    parseSEQs.pl --fasta - \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds.LRR-Pkinase.length.csv




##
## NB-ARC
##
awk '/NB-ARC/ {print $3;}' \
    /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.tbl | sort | uniq \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.NB-ARC.txt

fasta_process.pl --query /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.NB-ARC.txt \
    --fasta /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds.fa --rows 0 | \
    parseSEQs.pl --fasta - \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds.NB-ARC.length.csv


map_records.pl --query /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.LRR.txt \
    --subject /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.NB-ARC.txt \
    --rows1 0 --rows2 0 | awk '$2 != "N/A" {print $1;}' \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.NB-LRR.txt

fasta_process.pl --query /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.pfam.NB-LRR.txt \
    --fasta /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds.fa --rows 0 | \
    parseSEQs.pl --fasta - \
    > /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds.NB-LRR.length.csv
######################################################################




######################################################################
## gene features

cp -v /home/wl/Data/data/salix/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
    /home/wl/Data/biosoft/databases/snpEff/data/Ppersica_298_v2.1/sequences.fa
cp -v /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.cds.fa \
    /home/wl/Data/biosoft/databases/snpEff/data/Ppersica_298_v2.1/cds.fa
cp -v /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.transcript.fa \
    /home/wl/Data/biosoft/databases/snpEff/data/Ppersica_298_v2.1/gene.fa
cp -v /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.protein.fa \
    /home/wl/Data/biosoft/databases/snpEff/data/Ppersica_298_v2.1/protein.fa
cp -v /home/wl/Data/data/peach/ref/JGI_v2.1/annotation/Ppersica_298_v2.1.gene_exons.gff3 \
    /home/wl/Data/biosoft/databases/snpEff/data/Ppersica_298_v2.1/genes.gff

java -jar /home/wl/Data/biosoft/snpEff-4.0/snpEff.jar build \
    -config /home/wl/Data/biosoft/snpEff-4.0/snpEff.config \
    -gff3 -v Ppersica_298_v2.1 2>&1 | \
    tee /home/wl/Data/biosoft/databases/snpEff/data/Ppersica_298_v2.1/build.log


######################################################################


echo `date` ':all processes completed!'
echo

echo ---------------------------------------

exit 0




