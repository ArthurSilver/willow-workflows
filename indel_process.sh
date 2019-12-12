##generate initial candiate sites
bgzip -dc /home/ryf/Data/willow/ug/YAF1_s55.bwa.sort.dedup.realn.ug.vcf.gz | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/ins/ || /del/ || /\*/); my ($chrom, $pos, $ref)=(split /\s+/)[0,1,3];
        my $start = $pos-1; my $end = $pos+length($ref); print "$chrom\t$start\t$end\n";' | uniq \
    > /home/ryf/Data/willow/indel/ug/YAF1_s55.bwa.sort.dedup.realn.ug.indel.bed

bgzip -dc /home/ryf/Data/willow/hc/YAF1_s55.bwa.sort.dedup.realn.hc.vcf.gz | \
    vcf-annotate --fill-type | \
    perl -ne 'next if(/\#/); next unless(/ins/ || /del/ || /\*/); my ($chrom, $pos, $ref)=(split /\s+/)[0,1,3];
        my $start = $pos-1; my $end = $pos+length($ref); print "$chrom\t$start\t$end\n";' | uniq \
    > /home/ryf/Data/willow/indel/hc/YAF1_s55.bwa.sort.dedup.realn.hc.indel.bed

##recall use HC
Samples=`cat /home/ryf/Data/willow/indel/script/bamlist.txt | \
xargs -n 1 -P 4 -I PREFIX echo -n "-I PREFIX "`

/home/wl/bin/java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
-T HaplotypeCaller -nct 8 -stand_call_conf 30.0 -ip 10 \
-R /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
-L /home/ryf/Data/willow/indel/ug/YAF1_s55.bwa.sort.dedup.realn.ug.indel.bed \
-o /home/ryf/Data/willow/indel/ug/YAF1_s55.bwa.sort.dedup.realn.ug.indel.hc.vcf \
${Samples} 2>&1 | \
tee /home/ryf/Data/willow/indel/ug/YAF1_s55.bwa.sort.dedup.realn.ug.indel.hc.log

Samples=`cat /home/ryf/Data/willow/indel/script/bamlist.txt | \
xargs -n 1 -P 4 -I PREFIX echo -n "-I PREFIX "`

/home/wl/bin/java -jar /opt/nfs/share/biosoft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
-T HaplotypeCaller -nct 8 -stand_call_conf 30.0 -ip 10 \
-R /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.fa \
-L /home/ryf/Data/willow/indel/hc/YAF1_s55.bwa.sort.dedup.realn.hc.indel.bed \
-o /home/ryf/Data/willow/indel/hc/YAF1_s55.bwa.sort.dedup.realn.hc.indel.hc.vcf \
${Samples} 2>&1 | \
tee /home/ryf/Data/willow/indel/hc/YAF1_s55.bwa.sort.dedup.realn.hc.indel.hc.log

## screen process
##UG
### screen according to topology (by mutants)
vcf_process.pl --vcf /home/ryf/Data/willow/indel/ug/YAF1_s55.bwa.sort.dedup.realn.ug.indel.hc.vcf \
    --quality 30 --var-type indel | awk '/\#/ || ($9 ~ /AD/)' | \
    detect_mutations.pl -v - --max-cmp-total 3 --max-cmp-depth 2 --min-supp-depth 5 --max-cmp-miss 5 \
    -g /home/ryf/Data/willow/indel/script/groups.txt | \
    vcf-annotate -f c=2,150 --fill-type | \
    maskVCF.pl --input - --seq /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.softmasked.fa --add-pass \
    > /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.ug.indel.mut.c2t3d5m5.topology.vcf

### screen according to frequency
vcf_process.pl --vcf /home/ryf/Data/willow/indel/ug/YAF1_s55.bwa.sort.dedup.realn.ug.indel.hc.vcf \
    --quality 30 --var-type indel | awk '/\#/ || ($9 ~ /AD/)' | \
    detect_mutations.pl -v - --max-cmp-total 3 --max-cmp-depth 2 --min-supp-depth 5 --max-cmp-miss 0 --max-shared-freq 44 | \
    vcf-annotate -f c=2,150 --fill-type | \
    maskVCF.pl --input - --seq /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.softmasked.fa --add-pass \
    >  /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.ug.indel.mut.s44c2t3d5m0.frequency.vcf

##HC

### screen according to topology (by mutants)
vcf_process.pl --vcf /home/ryf/Data/willow/indel/hc/YAF1_s55.bwa.sort.dedup.realn.hc.indel.hc.vcf \
    --quality 30 --var-type indel | awk '/\#/ || ($9 ~ /AD/)' | \
    detect_mutations.pl -v - --max-cmp-total 3 --max-cmp-depth 2 --min-supp-depth 5 --max-cmp-miss 5 \
    -g /home/ryf/Data/willow/indel/script/groups.txt | \
    vcf-annotate -f c=2,150 --fill-type | \
    maskVCF.pl --input - --seq /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.softmasked.fa --add-pass \
    > /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.hc.indel.mut.c2t3d5m5.topology.vcf

### screen according to frequency
vcf_process.pl --vcf /home/ryf/Data/willow/indel/hc/YAF1_s55.bwa.sort.dedup.realn.hc.indel.hc.vcf \
    --quality 30 --var-type indel | awk '/\#/ || ($9 ~ /AD/)' | \
    detect_mutations.pl -v - --max-cmp-total 3 --max-cmp-depth 2 --min-supp-depth 5 --max-cmp-miss 0 --max-shared-freq 44 | \
    vcf-annotate -f c=2,150 --fill-type | \
    maskVCF.pl --input - --seq /home/hyw/Data/willow/ref/Spurpurea/assembly/Spurpurea_289_v1.0.softmasked.fa --add-pass \
    >  /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.hc.indel.mut.s44c2t3d5m0.frequency.vcf

##merge all results
## topology (by mutants)
vcf_process.pl --vcf /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.hc.indel.mut.c2t3d5m5.topology.vcf \
    --secondary-vcf /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.ug.indel.mut.c2t3d5m5.topology.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.combined.indel.mut.c2t3d5m5.topolgy.vcf

## frequency
vcf_process.pl --vcf /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.hc.indel.mut.s44c2t3d5m0.frequency.vcf \
    --secondary-vcf  /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.ug.indel.mut.s44c2t3d5m0.frequency.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag HC_GVCF --secondary-tag UG_Single --intersect-tag "UG_Single+HC_GVCF" \
    > /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.combined.indel.mut.s44c2t3d5m0.frequency.vcf

## combine topology and frequency
vcf_process.pl --vcf /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.combined.indel.mut.c2t3d5m5.topolgy.vcf \
    --secondary-vcf /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.combined.indel.mut.s44c2t3d5m0.frequency.vcf \
    --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag Grouped --secondary-tag NonGrouped --intersect-tag "Grouped+NonGrouped" \
    > /home/ryf/Data/willow/indel/combined/YAF1_s55.bwa.sort.dedup.realn.indel.mut.combined.vcf
