##UnifiedGenotyper call 
find /index/  -name "*.bwa.sort.dedup.realn.bam" -print | \
     xargs -n 1 -P 5 -I BAM_FILE sh -c '
        sample_id=`basename BAM_FILE | cut -d"." -f1`
        
        echo ">> Start calling ${sample_id} ... "
       
        java -jar GenomeAnalysisTK.jar \
          -R /ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
          -T UnifiedGenotyper -glm BOTH -nt 8 -stand_call_conf 30.0 -stand_emit_conf 30.0 -rf MappingQuality -mmq 20 \
          -o /index/ug/${sample_id}.bwa.sort.dedup.realn.ug.vcf \
          -I BAM_FILE \
          > /index/ug/${sample_id}.bwa.sort.dedup.realn.vcf.log 2>&1
     
        echo ">> Finished calling ${sample_id}"  
     '

ls /index/*.bwa.sort.dedup.realn.ug.vcf | \
    xargs -n 1 -P 4 -I {} /opt/nfs/bin/bgzip {}
ls /index/*.bwa.sort.dedup.realn.ug.vcf.gz | \
    xargs -n 1 -P 4 -I {} /opt/nfs/bin/tabix -p vcf {}

##HaplotypeCaller call
##variant discovery in GVCF mode
find /index/  -name "*.bwa.sort.dedup.realn.bam" -print | \
     xargs -n 1 -P 5 -I BAM_FILE sh -c '
        sample_id=`basename BAM_FILE | cut -d"." -f1`
        
        echo ">> Start calling ${sample_id} ... "
        
        java -jar GenomeAnalysisTK.jar \
            -R /ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
            -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
            -o /index/${sample_id}.bwa.sort.dedup.realn.hc.gvcf \
            -I BAM_FILE \
            > /index/${sample_id}.bwa.sort.dedup.realn.hc.log 2>&1
        
        echo ">> Finished calling ${sample_id}"
    '


ls/index/*.bwa.sort.dedup.realn.hc.gvcf | \
    xargs -n 1 -P 4 -I {} /opt/nfs/bin/bgzip {}
ls/index/*.bwa.sort.dedup.realn.hc.gvcf.gz | \
    xargs -n 1 -P 4 -I {} /opt/nfs/bin/tabix -p vcf {}
