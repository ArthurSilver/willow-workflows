######################################################################
## stats of depth and coverage


find /home/hyw/Data/willow/02.assembly/01.processed/ -name "*.bwa.sort.dedup.realn.bam" | \
    sed 's/.bam$//' | xargs -n 1 -P 1 -I PREFIX \
    sh -c '
        echo "[`date`]: Start processing PREFIX"
        
        samtools idxstats PREFIX.bam > PREFIX.idxstats
        
        echo "[`date`]: Finished processing PREFIX"
    '

grep "" /home/hyw/Data/willow/02.assembly/01.processed/*.bwa.sort.dedup.realn.idxstats \
    > /home/hyw/Data/willow/02.assembly/01.processed/YAF1_s48.bwa.sort.dedup.realn.idxstats.csv

## check GC contents
##
find /home/hyw/Data/willow/02.assembly/01.processed/ -name "*.bwa.sort.dedup.realn.bam" | \
    sed 's/.bam$//' | xargs -n 1 -P 4 -I PREFIX \
    sh -c '
        echo "[`date`]: Start processing PREFIX"
        
        java -jar /mnt/hp/wl/Data/biosoft/picard-tools-1.114/CollectGcBiasMetrics.jar \
            REFERENCE_SEQUENCE=/home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
            VALIDATION_STRINGENCY=LENIENT \
            INPUT=PREFIX.bam OUTPUT=PREFIX.GcBiasMetrics.csv \
            CHART_OUTPUT=PREFIX.GcBiasMetrics.pdf \
            SUMMARY_OUTPUT=PREFIX.GcBiasMetrics.summary.csv \
            > PREFIX.GcBiasMetrics.log 2>&1
        
        echo "[`date`]: Finished processing PREFIX"
    '
    
##
## depth of coverage
##
BAM_FILEs=`find /home/hyw/Data/willow/02.assembly/01.processed/ -name "*.bwa.sort.dedup.realn.bam" -print | \
    xargs -I BAM_FILE echo -n "-I BAM_FILE "`
java -jar /mnt/hp/wl/Data/biosoft/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
    -R /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
    -T DepthOfCoverage -nt 4 -rf BadCigar \
    -omitBaseOutput -omitIntervals -omitLocusTable --nBins 99 --start 1 --stop 100 -mmq 20 \
    -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 -ct 40 -ct 50 -ct 100 \
    -o /home/hyw/Data/willow/02.assembly/01.processed/YAF1_s48.bwa.sort.dedup.realn.MQ20.depth \
    ${BAM_FILEs} 2>&1 | \
    tee /home/hyw/Data/willow/02.assembly/01.processed/YAF1_s48.bwa.sort.dedup.realn.MQ20.depth.log