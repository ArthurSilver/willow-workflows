bcftools merge -O v -l /home/ryf/Data/willow/snp/script/uglist.txt | \
sed 's/\.\/\./0\/0/g' | bgzip -c > /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.vcf.gz && \
tabix -p vcf /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.vcf.gz

Samples=`cat /home/ryf/Data/willow/snp/script/hclist.txt | \
    xargs -I GVCF_FILE echo -n "-V GVCF_FILE "`

/home/wl/bin/java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
    -R /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
    -T GenotypeGVCFs -nt 1 -stand_call_conf 30.0 \
    -o /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.vcf \
    ${Col17_Samples} 2>& 1 | \
    tee /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.log

bgzip /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.vcf && \
    tabix -p vcf /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.vcf.gz
    
    
## Step1: Generate initial candidate targets
## UnifiedGenotyper
##
bgzip -dc /home/ryf/Data/willow/ug/YAF1_s55.bwa.sort.dedup.realn.ug.vcf.gz | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/snp/); my ($chrom, $pos)=(split /\s+/)[0,1]; print "$chrom\t",($pos-1),"\t$pos\n";' \
    > /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.bed
    
## HaplotypeCaller
##
bgzip -dc /home/ryf/Data/willow/hc/YAF1_s55.bwa.sort.dedup.realn.hc.vcf.gz | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/snp/); my ($chrom, $pos)=(split /\s+/)[0,1]; print "$chrom\t",($pos-1),"\t$pos\n";' \
    > /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.bed
  
## merge candidate target regions
##
cat /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.bed \
    /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.bed | \
    sort -k1,1 -k2,3n | bedtools merge -i - \
    >  /home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.bed
    
######################################################################
## Step2: Count accurate allele depths for each locus and each sample

##
## SNV sites
##

## Base Alignment Quality (BAQ) is a new concept deployed in samtools-0.1.9+.
## It aims to provide an efficient and effective way to rule out false SNPs caused by nearby INDELs.
## The default settings for mpileup is to filter reads with bitwise flag 0X704.
## So for pileup generation the following reads will not considered at all from the bam files:
##  1) 0x0400 (aka 1024 or "d") duplicate
##  2) 0x0200 (aka 512 or "f") failed QC
##  3) 0x0100 (aka 256 or "s") non primary alignment
##  4) 0x0004 (aka 4 or "u") unmapped
## Apply -A to use anomalous read pairs in mpileup, which are not used by default (requring r874+).

cat bamlist.txt | sed 's/.bam$//' | xargs -n 1 -P 8 -I PREFIX \
    sh -c '
        sample=`basename PREFIX | cut -d"." -f1`
        
        echo "[`date`]: Start processing ${sample} ... "
        
        ## count in anomalous reads, mapping quality >= 20
        samtools mpileup -Ad100000 -q20 -f /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -l /home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /home/ryf/Data/willow/snp/readcounts/ARMQ20/YAF1_s55.bwa.sort.dedup.realn.snp.${sample}.MQ20.AR.mpileup
        
        java -jar /opt/nfs/share/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
             /home/ryf/Data/willow/snp/readcounts/ARMQ20/YAF1_s55.bwa.sort.dedup.realn.snp.${sample}.MQ20.AR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /home/ryf/Data/willow/snp/readcounts/ARMQ20/YAF1_s55.bwa.sort.dedup.realn.snp.${sample}.MQ20.AR.readcounts
        
        ## count in anomalous reads, mapping quality >= 0
        samtools mpileup -Ad100000 -q0 -f /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -l /home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /home/ryf/Data/willow/snp/readcounts/ARMQ0/YAF1_s55.bwa.sort.dedup.realn.snp.${sample}.MQ0.AR.mpileup
        
        java -jar /opt/nfs/share/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
             /home/ryf/Data/willow/snp/readcounts/ARMQ0/YAF1_s55.bwa.sort.dedup.realn.snp.${sample}.MQ0.AR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /home/ryf/Data/willow/snp/readcounts/ARMQ0/YAF1_s55.bwa.sort.dedup.realn.snp.${sample}.MQ0.AR.readcounts
        
        ## only use proper pairs, mapping quality >= 20
        samtools mpileup -d100000 -q20 -f /opt/nfs/share/data/ref/arabidopsis/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -l /home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /home/ryf/Data/willow/snp/readcounts/NARMQ20/YAF1_s55.bwa.sort.dedup.realn.snp.${sample}.MQ20.NAR.mpileup
        
        java -jar /opt/nfs/share/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
            /home/ryf/Data/willow/snp/readcounts/NARMQ20/YAF1_s55.bwa.sort.dedup.realn.snp.${sample}.MQ20.NAR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /home/ryf/Data/willow/snp/readcounts/NARMQ20/YAF1_s55.bwa.sort.dedup.realn.snp.${sample}.MQ20.NAR.readcounts
        
        echo "[`date`]: Finished processing ${sample}"
    '


## count in anomalous reads, mapping quality >= 20
for f in `find /home/ryf/Data/willow/snp/readcounts \
    -name "YAF1_s55.bwa.sort.dedup.realn.snp.*.MQ20.AR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f8`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.MQ20.AR.readcounts.list

## count in anomalous reads, mapping quality >= 0
for f in `find /home/ryf/Data/willow/snp/readcounts \
    -name "YAF1_s55.bwa.sort.dedup.realn.snp.*.MQ0.AR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f8`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.MQ0.AR.readcounts.list

## count proper mapped reads only, mapping quality >= 20
for f in `find /\/home/ryf/Data/willow/snp/readcounts \
    -name "YAF1_s55.bwa.sort.dedup.realn.snp.*.MQ20.NAR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f8`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.MQ20.NAR.readcounts.list


##update AD
find /home/ryf/Data/willow/snp/combined/ \
    -name "YAF1_s55.bwa.sort.dedup.realn.snp.*.readcounts.list" | \
    xargs -n 1 -P 3 -I INPUT \
    sh -c '
        out_set=`echo INPUT | sed "s/.*.snp.//" | sed "s/.readcounts.list//"`
        
        echo "${out_set}"
        
        ## UnifiedGenotyper
        fillVcfDepth.pl --vcf /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.vcf.gz \
            --list INPUT --minimum-vcf --update-AD | bgzip -c \
            > /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.${out_set}.vcf.gz && \
            tabix -p vcf /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.${out_set}.vcf.gz
        
        ## HaplotypeCaller
        fillVcfDepth.pl --vcf /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.vcf.gz \
            --list INPUT --minimum-vcf --update-AD | bgzip -c \
            > /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.${out_set}.vcf.gz && \
            tabix -p vcf /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.${out_set}.vcf.gz
    '

######################################################################
## Step3: Screen out candidate mutations

## substitutions
find /home/ryf/Data/willow/ -name "YAF1_s55.bwa.sort.dedup.realn.*.snp.*.vcf.gz" | \
    xargs -n 1 -P 4 -I VCFIN \
    sh -c '
        out_prefix=`basename VCFIN | sed "s/vcf.gz/mut/"`
        out_prefix="/home/ryf/Data/willow/snp/${out_prefix}"
        
        echo "${out_prefix}"
        
        ### screen according to topology
        detect_mutations.pl -v VCFIN --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
            -g /home/ryf/Data/willow/snp/script/groups.txt | \
            vcf-annotate -f c=3,150 --fill-type | \
            maskVCF.pl --input - --seq /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.softmasked.fa --add-pass \
            --output ${out_prefix}.c2t3d5m5.topo.vcf
        
        ### screen according to frequency
        detect_mutations.pl -v VCFIN \
            --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 0 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 44 | \
            vcf-annotate -f c=3,150 --fill-type | \
            maskVCF.pl --input - --seq /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.softmasked.fa --add-pass \
            --output ${out_prefix}.s44c2t3d5m0.freq.vcf
     '

######################################################################
## Step4: Collect all cadidate mutations from different sources

##
## base substitution
##


## combine different readcount sets

### topology (by individuals)

#### UnifiedGenotyper
vcf_process.pl --vcf /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.mut.c2t3d5m5.topo.vcf \
    --secondary-vcf /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.MQ20.AR.mut.c2t3d5m5.topo.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.MQ0.AR.mut.c2t3d5m5.topo.vcf \
    > /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.mut.c2t3d5m5.topo.combined.vcf

#### HaplotypeCaller
vcf_process.pl --vcf /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.MQ20.NAR.mut.c2t3d5m5.topo.vcf \
    --secondary-vcf /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.MQ20.AR.mut.c2t3d5m5.topo.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.MQ0.AR.mut.c2t3d5m5.topo.vcf \
    > /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.mut.c2t3d5m5.topo.combined.vcf

vcf_process.pl --vcf /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.mut.c2t3d5m5.topo.combined.vcf \
    --secondary-vcf /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.mut.c2t3d5m5.topo.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    >/home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.mut.c2t3d5m5.topo.combined.vcf

### frequency

#### UnifiedGenotyper
vcf_process.pl --vcf /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.MQ20.NAR.mut.s44c2t3d5m0.freq.vcf \
    --secondary-vcf /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.MQ20.AR.mut.s44c2t3d5m0.freq.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.MQ0.AR.mut.s44c2t3d5m0.freq.vcf \
    > /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.mut.s44c2t3d5m0.freq.combined.vcf

#### HaplotypeCaller
vcf_process.pl --vcf /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.MQ20.NAR.mut.s44c2t3d5m0.freq.vcf \
    --secondary-vcf /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.MQ20.AR.mut.s44c2t3d5m0.freq.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.MQ0.AR.mut.s44c2t3d5m0.freq.vcf \
    > /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.mut.s44c2t3d5m0.freq.combined.vcf

vcf_process.pl --vcf /home/ryf/Data/willow/snp/hc/YAF1_s55.bwa.sort.dedup.realn.hc.snp.mut.s44c2t3d5m0.freq.combined.vcf \
    --secondary-vcf /home/ryf/Data/willow/snp/ug/YAF1_s55.bwa.sort.dedup.realn.ug.snp.mut.s44c2t3d5m0.freq.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    >/home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.mut.s44c2t3d5m0.freq.combined.vcf


## combine topology and frequency
vcf_process.pl --vcf /home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.mut.c2t3d5m5.topo.combined.vcf \
    --secondary-vcf /home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.mut.s44c2t3d5m0.freq.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag Grouped --secondary-tag NonGrouped --intersect-tag "Grouped+NonGrouped" \
    > /home/ryf/Data/willow/snp/combined/YAF1_s55.bwa.sort.dedup.realn.snp.mut.combined.vcf
#

