mkdir /path/readcounts
mkdir /path/combined
mkdir /path/hc_gvcf
mkdir /path/mapped
mkdir /path/reads
mkdir /path/ug_single
mkdir /path/variants

## simulate reads, reads in following situations were removed (FLAG 3844):
##    read unmapped
##    not primary alignment
##    read fails platform/vendor quality checks
##    read is PCR or optical duplicate
##    supplementary alignment
##
sim_mutation_reads.pl --fasta sequence.filtered.fasta \
    --depth /path/Plant_s21.bwa.sort.dedup.realn.mut.leaf.sampling_depths.csv \
    --random-size 1000 --samtools "-F 3844" --exclude mitochondria chloroplast \
    --group-file /path/script/groups-leaf.txt \
    > /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.dat

sim_mutation_reads.pl --fasta sequence.filtered.fasta \
    --depth /path/Plant_s21.bwa.sort.dedup.realn.mut.root.sampling_depths.csv \
    --random-size 1000 --samtools "-F 3844" --exclude mitochondria chloroplast \
    --group-file /path/script/groups-root.txt \
    >> /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.dat


## extract all reads cover the simulated regions in each sample
##
awk 'BEGIN{OFS="\t"} /^\#Tag/ || $1 == "MUT" {if(/\#Tag/){$2 = "#Chrom";} print;}' \
    /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.dat | \
    cut -f 2- | body sort -k1,1 -k2,2n \
    > /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.vars.csv

awk 'BEGIN{OFS="\t";} {if(/#Chrom/) {print "##fileformat=VCFv4.1";
        print "##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Simulated read depth\">";
        print "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Variant type\">";
        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUTATION";next;}
        print $1,$2,$4,$5,$6,"999\t.\tTYPE=snp\tSDP",$7;}' \
    /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.vars.csv \
    > /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.vars.vcf

awk '$1 ~ /SAM:/' /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.dat | \
    cut -f 2- \
    > /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.reads.sam

cat bamlist.txt | \
    xargs -n 1 -P 8 -I BAM_FILE sh -c '
        sample=`basename BAM_FILE | cut -d"." -f1`
        
        echo "${sample}"
        
        extract_bam_pairs.pl --samtools "-F 3844" --extend 2000 --bam BAM_FILE \
            --input /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.vars.csv \
            --patches /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.reads.sam \
            > /path/reads/${sample}_ex2k.sam \
            2> /path/reads/${sample}_ex2k.log
        
        samtools view -H BAM_FILE | \
            cat - /path/reads/${sample}_ex2k.sam \
            > /path/reads/${sample}_ex2k.add.sam
        
        ## sort bam file
        java -Djava.io.tmpdir=/home/ryf/Data/tmp -jar /opt/nfs/share/biosoft/picard-tools-1.114/SortSam.jar \
            INPUT=/path/reads/${sample}_ex2k.add.sam \
            OUTPUT=/path/reads/${sample}_ex2k.sort.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
            > /path/reads/${sample}_ex2k.sort.log 2>&1
        
        samtools index /path/reads/${sample}_ex2k.sort.bam
        
        sam2fastq.pl -i /path/reads/${sample}_ex2k.sam \
            -o /path/reads/${sample}_ex2k \
            >> /path/reads/${sample}_ex2k.log 2>&1 && \
            rm -v /path/reads/${sample}_ex2k.sam \
                  /path/reads/${sample}_ex2k.add.sam
    '

find /path/reads -name "*.fq" | \
    xargs -n 1 -P 12 -I {} gzip -v {}

## mapping and pre-processes

find /path/reads/ -name "*_1.fq.gz" | \
    sed 's/_1.fq.gz$//' | xargs -n 1 -P 8 -I PREFIX \
    sh -c '
        sample=`basename PREFIX | sed "s/_ex2k//"`
        lane_id=`basename PREFIX | sed "s/_ex2k//"`
        
        echo "[`date`]: Start mapping ${sample}:${lane_id} ... "
        
        read1=PREFIX"_1.fq.gz"
        read2=PREFIX"_2.fq.gz"
        
        ## Align reads with BWA-MEM algorithm
        bwa mem -t 1 -M -R "@RG\tID:${lane_id}\tLB:${sample}\tPL:Illumina\tPU:${sample}\tSM:${sample}" \
           sequence.filtered.fasta ${read1} ${read2} \
            > /path/mapped/${sample}.NEWREF.bwa.sam \
            2> /path/mapped/${sample}.NEWREF.bwa.process.log
        
        ## sort bam file
        java -Djava.io.tmpdir=/home/ryf/Data/tmp -jar /opt/nfs/share/biosoft/picard-tools-1.114/SortSam.jar \
            INPUT=/path/mapped/${sample}.NEWREF.bwa.sam \
            OUTPUT=/path/mapped/${sample}.NEWREF.bwa.sort.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
            >> /path/mapped/${sample}.NEWREF.bwa.process.log 2>&1 && \
            rm -v /path/mapped/${sample}.NEWREF.bwa.sam
        
        echo "[`date`]: Start marking duplicates ${sample} ... "
        
        ## mark duplicates
        java -jar /opt/nfs/share/biosoft/picard-tools-1.114/MarkDuplicates.jar \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT \
            INPUT=/path/mapped/${sample}.NEWREF.bwa.sort.bam \
            OUTPUT=/path/mapped/${sample}.NEWREF.bwa.sort.dedup.bam \
            METRICS_FILE=/path/mapped/${sample}.NEWREF.bwa.sort.dedup.metrics \
            >> /path/mapped/${sample}.NEWREF.bwa.process.log 2>&1 && \
            rm -v /path/mapped/${sample}.NEWREF.bwa.sort.bam
        
        ## index bam file
        samtools index /path/mapped/${sample}.NEWREF.bwa.sort.dedup.bam
        
        
        echo "[`date`]: Start realigning ${sample} ... "
        
        ## realignment
        java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R sequence.filtered.fasta \
            -T RealignerTargetCreator -nt 1 \
            -o /path/mapped/${sample}.NEWREF.bwa.sort.dedup.realn.intervals \
            -I /path/mapped/${sample}.NEWREF.bwa.sort.dedup.bam \
            >> /path/mapped/${sample}.NEWREF.bwa.process.log 2>&1
        
        java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R sequence.filtered.fasta \
            -T IndelRealigner \
            -targetIntervals /path/mapped/${sample}.NEWREF.bwa.sort.dedup.realn.intervals \
            -o /path/mapped/${sample}.NEWREF.bwa.sort.dedup.realn.bam \
            -I /path/mapped/${sample}.NEWREF.bwa.sort.dedup.bam \
            >> /path/mapped/${sample}.NEWREF.bwa.process.log 2>&1 && \
            rm -v /path/mapped/${sample}.NEWREF.bwa.sort.dedup.bam \
                  /path/mapped/${sample}.NEWREF.bwa.sort.dedup.bam.bai
        
        echo "[`date`]: Finished processing ${sample}"
    '

## UnifiedGenotyper
find /path/mapped/ -name "*.realn.bam" | \
    xargs -n 1 -P 8 -I BAM_FILE \
    sh -c '
        sample_id=`basename BAM_FILE | cut -d"." -f1`
        
        echo "[`date`]: Start calling ${sample_id} ... "
        
        java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R sequence.filtered.fasta \
            -T UnifiedGenotyper -glm BOTH -nt 4 -stand_call_conf 30.0 -rf MappingQuality -mmq 20 \
            -o /path/variants/${sample_id}.NEWREF.bwa.sort.dedup.realn.ug.vcf \
            -I BAM_FILE \
            > /path/variants/${sample_id}.NEWREF.bwa.sort.dedup.realn.ug.log 2>&1
        
        echo "[`date`]: Finished calling ${sample_id} ... "
    '

find /path/variants/ -name "*.vcf" -print | \
    xargs -n 1 -P 4 -I {} bgzip {}
find /path/variants/ -name "*.vcf.gz" -print | \
    xargs -n 1 -P 4 -I {} tabix -p vcf {}

## merge vcf files
bcftools merge -O v /path/variants/*.NEWREF.bwa.sort.dedup.realn.ug.vcf.gz | \
    sed 's/\.\/\./0\/0/g' \
    > /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.vcf



## HaplotypeCaller
find /path/mapped/ -name "*.realn.bam" -print | \
    xargs -n 1 -P 8 -I BAM_FILE sh -c '
        sample_id=`basename BAM_FILE | cut -d"." -f1`
        
        echo "[`date`]: Start calling ${sample_id} ... "
        
        java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
            -R sequence.filtered.fasta \
            -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
            -o /path/variants/${sample_id}.NEWREF.bwa.sort.dedup.realn.hc.gvcf \
            -I BAM_FILE \
            > /path/variants/${sample_id}.NEWREF.bwa.sort.dedup.realn.hc.log 2>&1
        
        echo "[`date`]: Finished calling ${sample_id} ... "
    '

## Joint genotyping
find /path/variants/ -name "*.hc.gvcf" -print | \
    xargs -n 1 -P 4 -I {} bgzip {}
find /path/variants/ -name "*.hc.gvcf.gz" -print | \
    xargs -n 1 -P 4 -I {} tabix -p vcf {}

GVCF_FILEs=`find /path/variants/ -name "*.hc.gvcf.gz" -print | \
    xargs -I GVCF_FILE echo -n "-V GVCF_FILE "`
java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
    -R sequence.filtered.fasta \
    -T GenotypeGVCFs -nt 4 -stand_call_conf 30.0 \
    -o /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.vcf \
    ${GVCF_FILEs} 2>& 1 | \
    tee /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.log


## check how many sites are callable
vcf_process.pl --vcf /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.vcf \
    --secondary-vcf /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.vcf \
    --combine-rows 0 1 --compare-rows 2 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.combined.vcf

map_records.pl --subject /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.combined.vcf \
    --query /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.vars.csv \
    --rows1 0 1 --rows2 0 1 | cut -f 3,9 --complement | perl -ne 's/(.*)\s+.*?Combine=(.*?)\s+/$1\t$2\t/; print;' \
    > /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.vars.called.csv

## Detect mutations


##
## Step1: Generate initial candidate targets
##

## UnifiedGenotyper
vcf_process.pl --vcf /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.vcf \
    --quality 50 --rare-only 54 \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.vcf

cat /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.vcf | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/snp/); my ($chrom, $pos)=(split /\s+/)[0,1]; print "$chrom\t",($pos-1),"\t$pos\n";' \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.bed


## HaplotypeCaller
vcf_process.pl --vcf /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.vcf \
    --quality 50 --rare-only 54 \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.vcf

cat /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.vcf | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/snp/); my ($chrom, $pos)=(split /\s+/)[0,1]; print "$chrom\t",($pos-1),"\t$pos\n";' \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.bed

## merge candidate target regions
cat /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.bed \
    /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.bed | \
    sort -k1,1 -k2,3n | bedtools merge -i - \
    > /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.bed
## count reads of all alleles in each strand
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
find /path/mapped/ -name "*.NEWREF.bwa.sort.dedup.realn.bam" -print | \
    sed 's/.bam$//' | xargs -n 1 -P 4 -I PREFIX \
    sh -c '
        sample=`basename PREFIX | cut -d"." -f1`
        
        echo "[`date`]: Start processing ${sample} ... "
        
        ## count in anomalous reads, mapping quality >= 20
        samtools mpileup -Ad100000 -q20 -f sequence.filtered.fasta \
            -l /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /path/readcounts/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.${sample}.MQ20.AR.mpileup
        
        java -jar /opt/nfs/share/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
                /path/readcounts/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.${sample}.MQ20.AR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /path/readcounts/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.${sample}.MQ20.AR.readcounts
        
        
        ## count in anomalous reads, mapping quality >= 0
        samtools mpileup -Ad100000 -q0 -f sequence.filtered.fasta \
            -l /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /path/readcounts/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.${sample}.MQ0.AR.mpileup
        
        java -jar /opt/nfs/share/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
            /path/readcounts/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.${sample}.MQ0.AR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /path/readcounts/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.${sample}.MQ0.AR.readcounts
        
        ## only use proper pairs, mapping quality >= 20
        samtools mpileup -d100000 -q20 -f sequence.filtered.fasta \
            -l /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /path/readcounts/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.${sample}.MQ20.NAR.mpileup
        
        java -jar /opt/nfs/share/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
                /path/readcounts/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.${sample}.MQ20.NAR.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /path/readcounts/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.${sample}.MQ20.NAR.readcounts
        
        echo "[`date`]: Finished processing ${sample}"
    '


## count in anomalous reads, MQ20
for f in `find /path/readcounts \
    -name "salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.*.MQ20.AR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f10`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.MQ20.AR.readcounts.list

## count in anomalous reads, MQ0
for f in `find /path/readcounts/ \
    -name "salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.*.MQ0.AR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f10`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.MQ0.AR.readcounts.list


## count proper mapped reads only, MQ20
for f in `find /path/readcounts/ \
    -name "salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.*.MQ20.NAR.readcounts" | sort`;
do
    library=`basename $f | cut -d"." -f10`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.MQ20.NAR.readcounts.list

## UnifiedGenotyper

### with anomalous reads, MQ20
fillVcfDepth.pl --vcf /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.vcf \
    --list /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.MQ20.AR.readcounts.list \
    --minimum-vcf --update-AD \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.AR.vcf

### with anomalous reads, MQ0
fillVcfDepth.pl --vcf /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.vcf \
    --list /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.MQ0.AR.readcounts.list \
    --minimum-vcf --update-AD \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ0.AR.vcf

### without anomalous reads, MQ20
fillVcfDepth.pl --vcf /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.vcf \
    --list /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.MQ20.NAR.readcounts.list \
    --minimum-vcf --update-AD \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.NAR.vcf


## HaplotypeCaller

### with anomalous reads, MQ20
fillVcfDepth.pl --vcf /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.vcf \
    --list /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.MQ20.AR.readcounts.list \
    --minimum-vcf --update-AD \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.AR.vcf

### with anomalous reads, MQ0
fillVcfDepth.pl --vcf /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.vcf \
    --list /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.MQ0.AR.readcounts.list \
    --minimum-vcf --update-AD \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ0.AR.vcf

### without anomalous reads, MQ20
fillVcfDepth.pl --vcf /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.vcf \
    --list /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.MQ20.NAR.readcounts.list \
    --minimum-vcf --update-AD \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.NAR.vcf

## UnifiedGenotyper
##

## substitutions

### screen according to frequency

#### with anomalous reads, MQ>=20
detect_mutations.pl -v /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 44 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.AR.mut.s44c2t3d5m5.frequency.vcf

#### with anomalous reads, MQ>=0
detect_mutations.pl -v /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ0.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 44 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ0.AR.mut.s44c2t3d5m5.frequency.vcf

#### without anomalous reads, MQ>=20
detect_mutations.pl -v /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.NAR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 44 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.NAR.mut.s44c2t3d5m5.frequency.vcf

### screen according to topology
#### with anomalous reads, MQ>=20
detect_mutations.pl -v /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.AR.vcf \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 \
    -g /path/script/sample-group.txt | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.AR.mut.c2t3d5m5.topo.vcf

#### with anomalous reads, MQ>=0
detect_mutations.pl -v /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ0.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    -g /path/script/sample-group.txt \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ0.AR.mut.c2t3d5m5.topo.vcf

#### without anomalous reads, MQ>=20
detect_mutations.pl -v /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.NAR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    -g /path/script/sample-group.txt \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.NAR.mut.c2t3d5m5.topo.vcf

##
## HaplotypeCaller
##

## substitutions


### screen according to topology
#### with anomalous reads, MQ>=20
detect_mutations.pl -v /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    -g /path/script/sample-group.txt \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.AR.mut.c2t3d5m5.topo.vcf

#### with anomalous reads, MQ>=0
detect_mutations.pl -v /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ0.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    -g /path/script/sample-group.txt \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ0.AR.mut.c2t3d5m5.topo.vcf

#### without anomalous reads, MQ>=20
detect_mutations.pl -v /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.NAR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    -g /path/script/sample-group.txt \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.NAR.mut.c2t3d5m5.topo.vcf

### screen according to frequency

#### with anomalous reads, MQ>=20
detect_mutations.pl -v /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 44 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.AR.mut.s44c2t3d5m5.frequency.vcf

#### with anomalous reads, MQ>=0
detect_mutations.pl -v /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ0.AR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 44 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ0.AR.mut.s44c2t3d5m5.frequency.vcf

#### without anomalous reads, MQ>=20
detect_mutations.pl -v /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.NAR.vcf \
    --mask-only LowDepth StrandBias HighMissing \
    --max-cmp-depth 2 --max-cmp-total 3 --min-supp-depth 5 --max-cmp-miss 5 --min-supp-plus 1 --min-supp-minus 1 --max-shared-freq 44 | \
    vcf-annotate -f c=3,150 --fill-type | \
    maskVCF.pl --input - --seq /mnt/san2/usr/ryf/ref/Salix_suchowensis_v4.1/fasta/Ssuchowensis_4.fa.masked --add-pass \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.NAR.mut.s44c2t3d5m5.frequency.vcf

## Step4: Collect all cadidate mutations from different sources

##
## base substitution
##

## combine different readcount sets

### frequency

#### UnifiedGenotyper
vcf_process.pl --vcf /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.NAR.mut.s44c2t3d5m5.frequency.vcf \
    --secondary-vcf /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.AR.mut.s44c2t3d5m5.frequency.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ0.AR.mut.s44c2t3d5m5.frequency.vcf \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.mut.s44c2t3d5m5.frequency.combined.vcf

#### HaplotypeCaller
vcf_process.pl --vcf /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.NAR.mut.s44c2t3d5m5.frequency.vcf \
    --secondary-vcf /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.AR.mut.s44c2t3d5m5.frequency.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ0.AR.mut.s44c2t3d5m5.frequency.vcf \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.mut.s44c2t3d5m5.frequency.combined.vcf


## combine different methods

## frequency
vcf_process.pl --vcf /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.mut.s44c2t3d5m5.frequency.combined.vcf \
    --secondary-vcf /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.mut.s44c2t3d5m5.frequency.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.mut.s44c2t3d5m5.frequency.combined.vcf

#####topology

#### UnifiedGenotyper
vcf_process.pl --vcf /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.NAR.mut.c2t3d5m5.topo.vcf \
    --secondary-vcf /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ20.AR.mut.c2t3d5m5.topo.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.MQ0.AR.mut.c2t3d5m5.topo.vcf \
    > /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.mut.c2t3d5m5.topo.combined.vcf

#### HaplotypeCaller
vcf_process.pl --vcf /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.NAR.mut.c2t3d5m5.topo.vcf \
    --secondary-vcf /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ20.AR.mut.c2t3d5m5.topo.vcf \
    --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
    vcf_process.pl --vcf - --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
    --secondary-vcf /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.MQ0.AR.mut.c2t3d5m5.topo.vcf \
    > /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.mut.c2t3d5m5.topo.combined.vcf

## combine different methods
vcf_process.pl --vcf /path/hc_gvcf/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.hc.fq54.snp.mut.c2t3d5m5.topo.combined.vcf \
    --secondary-vcf /path/ug_single/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.ug.fq54.snp.mut.c2t3d5m5.topo.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.mut.c2t3d5m5.topo.combined.vcf

## combine topology and frequency
vcf_process.pl --vcf /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.mut.c2t3d5m5.topo.combined.vcf \
    --secondary-vcf /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.mut.s44c2t3d5m5.frequency.combined.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag Grouped --secondary-tag NonGrouped --intersect-tag "Grouped+NonGrouped" \
    > /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.snp.mut.combined.vcf
######################################################################

## analysis of results

## compare results
vcf_process.pl --vcf /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.snp.mut.combined.vcf \
    --secondary-vcf /path/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.simulated.vars.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag FPC --secondary-tag FNC --intersect-tag TPC \
    > /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.mut.vcf


perl -ne 'next if (/^\#\#/); if (/\#/) {
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUTATION\tMethods\tMut_Type\tFrequency\tMut_Set\tCategory\n"; next;}
    my @line = (split /\s+/); my $info = $line[7]; $info =~ /MA=(\w+)/; my $mut_allele = $1;
    my $fq_sum = 1; if($info =~ /Shared\=(\d+)/){$fq_sum = $1;}
    my $type = ($mut_allele eq $line[3]) ? "REF" : "ALT"; my $out_line = join "\t", @line;
    my $method = "UG_Single"; if ($info =~ /UG_Single\+HC_GVCF/){$method = "UG_Single+HC_GVCF";}elsif($info =~ /HC_GVCF/){$method = "HC_GVCF";}
    my @filters = (); if ($info =~ /Combine=AR;/) {push @filters, "AR";} if ($info =~ /Combine=NAR;/) {push @filters, "NAR";}
    if ($info =~ /Combine=MQ20;/) {push @filters, "MQ20";} if ($info =~ /Combine=MQ0;/) {push @filters, "MQ0";}
    if ($info !~ /NMISS=0/) {push @filters, "MISSING";} if (($info !~ /FPD=0/) && ($info !~ /FPD=1;FPFQ=1;/)) {push @filters, "FPD";}
    if ($method eq "UG_Single") {push @filters, "UG";} if ($line[5] < 50) {push @filters, "LowQual";}
    my $set = "TP"; if ($info =~ /Combine=Grouped\+NonGrouped/){$set .= "+FQ";}elsif($info =~ /Combine=NonGrouped/){$set = "FQ";}
    if (scalar @filters == 0) {$set .= "(Confidence)";} else {my $filters = join ",", @filters; $set .= "($filters)";}
    my $cate = "TPC"; if ($info =~ /FPC/){$cate = "FPC";}elsif ($info =~ /FNC/){$cate = "FNC";}
    print "$out_line\t$method\t$type\t$fq_sum\t$set\t$cate\n";' \
    /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.mut.vcf \
    > /path/combined/salix_s55.simu6.NEWREF.bwa.sort.dedup.realn.fq54.snp.mut.csv
