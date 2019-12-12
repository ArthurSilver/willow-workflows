##UnifiedGenotyper call 
find /home/hyw/Data/willow/02.assembly/process1/  -name "*.bwa.sort.dedup.realn.bam" -print | \
     xargs -n 1 -P 5 -I BAM_FILE sh -c '
        sample_id=`basename BAM_FILE | cut -d"." -f1`
        
        echo ">> Start calling ${sample_id} ... "
       
        java -jar /home/hyw/Data/software/GATK/GenomeAnalysisTK.jar \
          -R /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
          -T UnifiedGenotyper -glm BOTH -nt 8 -stand_call_conf 30.0 -stand_emit_conf 30.0 -rf MappingQuality -mmq 20 \
          -o /home/hyw/Data/willow/03.analysis/01.variants/01.multi/ug/${sample_id}.bwa.sort.dedup.realn.ug.vcf \
          -I BAM_FILE \
          > /home/hyw/Data/willow/03.analysis/01.variants/01.multi/ug/${sample_id}.bwa.sort.dedup.realn.vcf.log 2>&1
     
        echo ">> Finished calling ${sample_id}"  
     '

ls /home/hyw/Data/willow/03.analysis/01.variants/00.single/*.bwa.sort.dedup.realn.ug.vcf | \
    xargs -n 1 -P 4 -I {} /opt/nfs/bin/bgzip {}
ls /home/hyw/Data/willow/03.analysis/01.variants/00.single/*.bwa.sort.dedup.realn.ug.vcf.gz | \
    xargs -n 1 -P 4 -I {} /opt/nfs/bin/tabix -p vcf {}

##HaplotypeCaller call
##variant discovery in GVCF mode
find /home/hyw/Data/willow/02.assembly/process1/  -name "*.bwa.sort.dedup.realn.bam" -print | \
     xargs -n 1 -P 5 -I BAM_FILE sh -c '
        sample_id=`basename BAM_FILE | cut -d"." -f1`
        
        echo ">> Start calling ${sample_id} ... "
        
        java -jar /home/hyw/Data/software/GATK/GenomeAnalysisTK.jar \
            -R /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
            -et NO_ET -K /home/hyw/Data/software/evolution_smail.nju.edu.cn.key \
            -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
            -o /home/hyw/Data/willow/03.analysis/01.variants/00.single/${sample_id}.bwa.sort.dedup.realn.hc.gvcf \
            -I BAM_FILE \
            > /home/hyw/Data/willow/03.analysis/01.variants/00.single/${sample_id}.bwa.sort.dedup.realn.hc.log 2>&1
        
        echo ">> Finished calling ${sample_id}"
    '

find /home/hyw/Data/willow/02.assembly/process2/  -name "*.bwa.sort.dedup.realn.bam" -print | \
     xargs -n 1 -P 10 -I BAM_FILE sh -c '
        sample_id=`basename BAM_FILE | cut -d"." -f1`
        
        echo ">> Start calling ${sample_id} ... "
        
        java -jar /home/hyw/Data/software/GATK/GenomeAnalysisTK.jar \
            -R /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
            -et NO_ET -K /home/hyw/Data/software/evolution_smail.nju.edu.cn.key \
            -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
            -o /mnt/ibm4/hyw/willow//03.analysis/01.variants/00.single/${sample_id}.bwa.sort.dedup.realn.hc.gvcf \
            -I BAM_FILE \
            > /mnt/ibm4/hyw/willow/03.analysis/01.variants/00.single/${sample_id}.IBM5.bwa.sort.dedup.realn.hc.log 2>&1
        
        echo ">> Finished calling ${sample_id}"
    '


ls /home/hyw/Data/willow/03.analysis/01.variants/00.single/*.bwa.sort.dedup.realn.hc.gvcf | \
    xargs -n 1 -P 4 -I {} /opt/nfs/bin/bgzip {}
ls /home/hyw/Data/willow/03.analysis/01.variants/00.single/*.bwa.sort.dedup.realn.hc.gvcf.gz | \
    xargs -n 1 -P 4 -I {} /opt/nfs/bin/tabix -p vcf {}
