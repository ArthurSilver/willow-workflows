find -L /index/01.cleandata/clean_read2/ -name "*_1.fq.gz" | sed 's/_1.fq.gz$//' | \
    xargs -n 1 -P 2 -I PREFIX \
    sh -c '
        lane_id=`basename PREFIX`
        sample=`dirname PREFIX`
        sample=`basename ${sample}`
        
        echo "[`date`]: Start mapping ${sample}:${lane_id} ... "
        
        read1=PREFIX"_1.fq.gz"
        read2=PREFIX"_2.fq.gz"
        
        ## Align reads with BWA-MEM algorithm
       /bwa/bwa-0.7.12/bwa mem -t 3 -M -R "@RG\tID:${lane_id}\tLB:${sample}\tPL:Illumina\tPU:${sample}\tSM:${sample}" \
            /index/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa ${read1} ${read2} \
            > /index/02.assembly/00.mapped/${sample}.bwa.sam \
            2> /index/02.assembly/00.mapped/${sample}.bwa.log
        
        ## sort bam file
        java -Djava.io.tmpdir=/Data/tmp -jar /picard/picard-tools-1.124/picard.jar SortSam \
            INPUT=/index/02.assembly/00.mapped/${sample}.bwa.sam \
            OUTPUT=/index/02.assembly/00.mapped/${sample}.bwa.sort.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
            >> /index/02.assembly/00.mapped/${sample}.bwa.log 2>&1 && \
            rm -v /index/02.assembly/00.mapped/${sample}.bwa.sam
        
        echo "[`date`]: Start marking duplicates ${sample} ... "
        
        ## mark duplicates
        java -jar /picard/picard-tools-1.124/picard.jar MarkDuplicates \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT \
            INPUT=/index/02.assembly/00.mapped/${sample}.bwa.sort.bam \
            OUTPUT=/index/02.assembly/00.mapped/${sample}.bwa.sort.dedup.bam \
            METRICS_FILE=/index/02.assembly/00.mapped/${sample}.bwa.sort.dedup.metrics \
            >> /index/02.assembly/00.mapped/${sample}.bwa.log 2>&1 && \
            rm -v /index/02.assembly/00.mapped/${sample}.bwa.sort.bam
        
        ## index bam file
        samtools index /index/02.assembly/00.mapped/${sample}.bwa.sort.dedup.bam
        
        echo "[`date`]: Start realigning ${sample} ... "
        
        ## realignment
        java -jar /GATK/GenomeAnalysisTK.jar \
            -R /index/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
            -T RealignerTargetCreator -nt 2 \
            -o /index/02.assembly/00.mapped/${sample}.bwa.sort.dedup.realn.intervals \
            -I /index/02.assembly/00.mapped/${sample}.bwa.sort.dedup.bam \
            >> /index/02.assembly/00.mapped/${sample}.bwa.log 2>&1
        
        java -jar /GATK/GenomeAnalysisTK.jar \
            -R /index/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
            -T IndelRealigner \
            -targetIntervals /index/02.assembly/00.mapped/${sample}.bwa.sort.dedup.realn.intervals \
            -o /index/02.assembly/00.mapped/${sample}.bwa.sort.dedup.realn.bam \
            -I /index/02.assembly/00.mapped/${sample}.bwa.sort.dedup.bam \
            >> /index/02.assembly/00.mapped/${sample}.bwa.log 2>&1 && \
            rm -v /index/02.assembly/00.mapped/${sample}.bwa.sort.dedup.bam \
                  /index/02.assembly/00.mapped/${sample}.bwa.sort.dedup.bam.bai
        
        echo "[`date`]: Finished processing ${sample}"
    '