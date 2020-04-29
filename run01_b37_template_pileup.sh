##!/bin/bash
#
##USAGE: work.sh sample_name
#
sample=$1
threads="1"
main="/nfs/home/Mandy"
maindir="${main}/BJCH_redo_190510"
#fqpath="${main}/BJCH/BJCH_fq_original"
bamdir="${maindir}/BAM"
resultdir="${maindir}/new_result"
logdir="${maindir}/log"
#if [ ! -d ${outdir}/${sample} ];then
#	mkdir ${outdir}/${sample}
#fi
#resultdir=${outdir}/${sample}
ref37="${main}/resources/b37/ref_Homo_b37/human_g1k_v37.fasta"
ref37_2bit="${main}/resources/b37/ref_Homo_b37/human_g1k_v37.2bit"
dbsnp="${main}/resources/b37/ref_Homo_b37/dbsnp_138.b37.vcf"
mills="${main}/resources/b37/ref_Homo_b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
G1000snps="${main}/resources/b37/ref_Homo_b37/1000G_phase1.snps.high_confidence.b37.vcf"
G1000indels="${main}/resources/b37/ref_Homo_b37/1000G_phase1.indels.b37.vcf"
bwa="${main}/software/bwa/bwa"
picard="${main}/software/picard.jar"
gatk="${main}/software/gatk-4.0.4.0/gatk"
samtools="${main}/software/samtools-1.8/samtools"
annovar_path="${main}/software/annovar"
intervar_path="${main}/software/InterVar"
fastp="${main}/software/fastp"
fastqc="${main}/software/FastQC/./fastqc"
#ad_dp_filter_dp30="${main}/BJCH/bin/germline_ad_dp30_filter.py"
ad_dp_filter_dp30="${main}/BJCH/bin/germline_ad_dp_filter.py"
ad_dp_filter_dp10="${main}/BJCH/bin/germline_af0.2_dp10_filter.py"
af_dp_filter_split="${main}/BJCH/bin/germline_eachalt_af0.2dp10_filter.py"
vcf_split="${main}/BJCH/bin/vcf_split.py"
HC_SB_filter="${main}/BJCH/bin/germline_HC_SB_filter.py"
pileup_sum="${main}/BJCH/bin/mpileup_strandbias.py"
ad_dp_filter_dp20="${main}/BJCH/bin/software/germline_ad_dp20_filter.py"
#filter_intervar_LOF="${WES_main}/BC/bin/filter_intervar.py"
mpile_sum="${main}/BJCH/bin/mpile_sum.py"
#gnomad2_1_annotate_oea="${WES_main}/BC/bin/gnomad2.1_mysql_OEAS_LOF_annotate.py"
#gnomad2_0_AF="${WES_main}/BC/bin/gnomad_EAS_AF_0.001.py"
#gnomad2_1_AF="${WES_main}/BC/bin/gnomad2.1_EAS_AF_0.001.py"
depth_sum="${main}/BJCH/bin/depth_sum.pl"
prep_bed="${main}/BJCH/bin/prep_bed.pl"
#vcf_processor="${WES_main}/BC/bin/vcf_processor.py"
#check_zygo="${WES_main}/BC/bin/check_var_hom_het.py"
echo `date`
#echo "Now processing $sample"
#
###${cutadapt} -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
##	--error-rate=0.1 --trim-n \
##	--minimum-length=20 \
##	-o ${resultdir}/${sample}.R1.cutadapt3.fq.gz \
##	-p ${resultdir}/${sample}.R2.cutadapt3.fq.gz \
##	${fqpath}/${sample}.R1.clean.fastq.gz ${fqpath}/${sample}.R2.clean.fastq.gz > ${resultdir}/${sample}.cut3.summary
###
###${cutadapt} -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG \
##	--error-rate=0.1 --minimum-length=20 \
##	--trim-n --minimum-length=20 \
##	-o ${resultdir}/${sample}.R1.cutadapt35.fq.gz \
##	-p ${resultdir}/${sample}.R2.cutadapt35.fq.gz \
##	${resultdir}/${sample}.R1.cutadapt3.fq.gz ${resultdir}/${sample}.R2.cutadapt3.fq.gz > ${resultdir}/${sample}.cutadapt_35.summary
###
###${sickle} pe \
##	-f ${resultdir}/${sample}.R1.cutadapt35.fq.gz \
##	-r ${resultdir}/${sample}.R2.cutadapt35.fq.gz \
##	-t sanger -g \
##	-o ${resultdir}/${sample}.R1.trimmed.fq.gz \
##	-p ${resultdir}/${sample}.R2.trimmed.fq.gz \
##	-s ${resultdir}/${sample}.trimmed.fq.gz > ${resultdir}/${sample}.trim.summary
################FASTQC-ASSESSMENT###############
#${fastqc} -t ${threads} -o ${resultdir} \
#	${fqpath}/${sample}*R1.clean.fastq.gz \
#	${fqpath}/${sample}*R2.clean.fastq.gz > ${resultdir}/${sample}.fastqc.summary
#
#######FASTP-trimming#######
#${fastp} -w ${threads} -i ${fqpath}/${sample}*R1.clean.fastq.gz -I ${fqpath}/${sample}*R2.clean.fastq.gz \
#-o ${resultdir}/${sample}.R1.trimmed.fq.gz -O ${resultdir}/${sample}.R2.trimmed.fq.gz \
#--correction \
#--json ${resultdir}/${sample}.json \
#--html ${resultdir}/${sample}.html \
#--report_title ${sample}" fastp report" > ${resultdir}/${sample}.fastp.summary.txt &&
#
##########BWA-Alignment#########
#${bwa} mem \
#	-M -V -t ${threads} \
#	-R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina\tLB:${sample}\tPU:HWIorHISEQ" \
#	${ref37} \
#	${resultdir}/${sample}.R1.trimmed.fq.gz ${resultdir}/${sample}.R2.trimmed.fq.gz \
#	-o ${resultdir}/${sample}.aln.bam &&
#
##########PICARD#########
#java -jar -Xmx4g ${picard} SortSam \
#	I=${resultdir}/${sample}.aln.bam \
#	O=${resultdir}/${sample}.sort.bam \
#	SORT_ORDER=coordinate 
#
#java -jar -Xmx4g ${picard} MarkDuplicates \
#	I=${resultdir}/${sample}.sort.bam \
#	O=${resultdir}/${sample}.dedup.bam \
#	M=${resultdir}/${sample}.dedup_metrix.txt \
#	ASSUME_SORTED=true \
#	REMOVE_DUPLICATES=true \
#	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100
#
#java -jar -Xmx4g ${picard}  FixMateInformation \
#	I=${resultdir}/${sample}.dedup.bam \
#	O=${resultdir}/${sample}.fixmat.bam \
#	SO=coordinate AS=true
#
###########GATK##################
##${gatk} --java-options "Xmx10g" BQSRPipelineSpark \
##	-R ${ref37_2bit} \
##	-I ${resultdir}/${sample}.dedup.bam \
##	--known-sites ${dbsnp} \
##	--known-sites ${mills} \
##	--known-sites ${G1000snps} \
##	--know-sites ${G1000indels} \
##	--bqsr-baq-gap-open-penalty 30 \
##	--gcs-max-retries 100 \
##	-O ${resultdir}/${sample}.recal.bam &&
###################################
#${gatk} --java-options "-Xmx4g" BaseRecalibrator \
#	-R ${ref37} \
#	-I ${resultdir}/${sample}.fixmat.bam \
#	--known-sites ${dbsnp} \
#	--known-sites ${mills} \
#	--known-sites ${G1000snps} \
#	--known-sites ${G1000indels} \
#	-O ${resultdir}/${sample}.recal_date.grp &&
#
#${gatk} --java-options "-Xmx4g" ApplyBQSR \
#	-R ${ref37} \
#	-I ${resultdir}/${sample}.fixmat.bam \
#	-O ${resultdir}/${sample}.recal.bam \
#	--bqsr-recal-file ${resultdir}/${sample}.recal_date.grp &&
#
#
#${gatk} BuildBamIndex -I=${resultdir}/${sample}.recal.bam &&
#
#${gatk} --java-options "-Xmx4g" HaplotypeCaller \
#	--native-pair-hmm-threads ${threads} \
#	-R ${ref37} \
#	-I ${resultdir}/${sample}.recal.bam \
#	--dbsnp ${dbsnp} \
#  --annotation StrandBiasBySample \
#	-O ${resultdir}/${sample}.sb.raw.snps.indels.vcf \
#	-bamout ${resultdir}/${sample}.sb.final.bam &&
#
#${samtools} depth \
#	--reference ${ref37} \
#	${resultdir}/${sample}.sb.final.bam > ${resultdir}/${sample}.sb.final.depth
#
#perl ${depth_sum} ${resultdir}/${sample}.sb.final.depth ${resultdir}/${sample}.sb.final.depth.sum
#perl ${prep_bed} ${resultdir}/${sample}.sb.final.depth.sum 100 20 >  ${resultdir}/${sample}.sb.final.depth.100.20
#perl ${prep_bed} ${resultdir}/${sample}.sb.final.depth.sum 200 30 >  ${resultdir}/${sample}.sb.final.depth.200.30

#cat ${resultdir}/${sample}.sb.final.depth.100.20 | perl -e 'use POSIX;$len=0;$dep=0;while(<>){chomp;@list=split(/\t/,$_);$len+=$list[3];$dep+=($list[3]*$list[4]);}print ceil($dep/$len)."\n"' > ${resultdir}/${sample}.coverage.100len.20dep
#cat ${resultdir}/${sample}.sb.final.depth.200.30 | cut -f 5| perl -e '$a=0;$b=0;while(<>){chomp;$a+=$_;$b+=1;}print int($a/$b)."\n"' > ${resultdir}/${sample}.coverage1.200.30
#cat ${resultdir}/${sample}.sb.final.depth.200.30 | perl -e '$a=0;$b=0;while(<>){chomp;@c=split;$a+=$c[3];$b+=($c[3]*$c[4]);}print int($b/$a)."\n"' > ${resultdir}/${sample}.coverage2.200.30
 
##
#######Select SNPs#########
#${gatk} SelectVariants \
#	-R ${ref37} \
#	-V ${resultdir}/${sample}.sb.raw.snps.indels.vcf \
#	-select-type SNP \
#	-O ${resultdir}/${sample}.sb.raw_snps.vcf &&
#
#${gatk} VariantFiltration \
#	-R ${ref37} \
#	-V ${resultdir}/${sample}.sb.raw_snps.vcf \
#	--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 ||MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
#	--filter-name "my_snp_filter" -O ${resultdir}/${sample}.sb.tagged_snps.vcf &&
#
#${gatk} SelectVariants \
#	-R ${ref37} \
#	-V ${resultdir}/${sample}.sb.tagged_snps.vcf \
#	--exclude-filtered true \
#	-O ${resultdir}/${sample}.sb.filtered_snps.vcf
#
##python ${ad_dp_filter_dp20} ${resultdir}/${sample}.filtered_snps.vcf ${resultdir}/${sample}.filtered_addp20_snps.vcf
###########Select Indels###############
#${gatk} SelectVariants \
#	-R ${ref37} \
#	-V ${resultdir}/${sample}.sb.raw.snps.indels.vcf \
#	-select-type INDEL \
#	-O ${resultdir}/${sample}.sb.raw_indels.vcf &&
#
#${gatk} VariantFiltration \
#	-R ${ref37} \
#	-V ${resultdir}/${sample}.sb.raw_indels.vcf \
#	--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
#	--filter-name "my_indel_filter" \
#	-O ${resultdir}/${sample}.sb.tagged_indels.vcf &&
#
#${gatk} SelectVariants \
#	-R ${ref37} \
#	-V ${resultdir}/${sample}.sb.tagged_indels.vcf \
#	--exclude-filtered true \
#	-O ${resultdir}/${sample}.sb.filtered_indels.vcf
##
#
#java -jar ${picard} MergeVcfs \
#	I=${vcfdir}/${sample}.filtered_snps.vcf \
#	I=${vcfdir}/${sample}.filtered_indels.vcf \
#	O=${resultdir}/${sample}.filtered.merged.vcf
##
#python ${vcf_processor} \
#	${resultdir}/${sample}.filtered.merged.vcf \
#	${resultdir}/${sample}.filtered.merged.split.vcf \
#	${resultdir}/${sample}.filtered.merged.split > ${logdir}/${sample}.processed.log
##	${resultdir}/${sample}.filtered.merged.split > ${logdir}/${sample}.processed.log
##	
##grep -i "#" ${resultdir}/${sample}.sb.filtered_indels.vcf > ${resultdir}/${sample}.sb.head
##
##grep -i -v "#" ${resultdir}/${sample}.sb.filtered_indels.vcf > ${resultdir}/${sample}.sb.filtered_indels2.vcf
##
##grep -i -v "#" ${resultdir}/${sample}.sb.filtered_snps.vcf > ${resultdir}/${sample}.sb.filtered_snps2.vcf
##cat ${resultdir}/${sample}.sb.head ${resultdir}/${sample}.sb.filtered_snps2.vcf ${resultdir}/${sample}.sb.filtered_indels2.vcf > ${resultdir}/${sample}.sb.filtered_variants.vcf
##
##python ${vcf_split} ${resultdir}/${sample}.sb.filtered_variants.vcf ${resultdir}/${sample}.sb.split.vcf
##
###python ${ad_dp_filter_dp10} ${resultdir}/${sample}.sb.filtered_snps.vcf ${resultdir}/${sample}.sb.filtered_af0.2dp10_snps.vcf
##
###python ${ad_dp_filter_dp10} ${resultdir}/${sample}.sb.filtered_indels.vcf ${resultdir}/${sample}.sb.filtered_af0.2dp10_indels.vcf
##python ${af_dp_filter_split} ${resultdir}/${sample}.sb.split.vcf ${resultdir}/${sample}.sb.filtered.split_af0.2dp10.vcf
###
##python ${HC_SB_filter} ${resultdir}/${sample}.sb.filtered.split_af0.2dp10.vcf ${resultdir}/${sample}.sb.filtered.split_af0.2dp10_sb.vcf
######USE ANNOVAR to annotate the AF_EAS in gnomAD2.0###
##perl ${annovar_path}/table_annovar.pl \
##	${resultdir}/${sample}.filtered.merged.split.filterDPADAF.vcf \
##	${annovar_path}/humandb/ \
##	-buildver hg19 \
##	-out ${resultdir}/${sample}.processed_gnomAD \
##	-remove \
##	-protocol gnomad_exome \
##	-operation f \
##	-nastring . \
##	-vcfinput
#######USE gnomAD2.1 in mysql built in-house to annotate the AF_oea###
#perl ${annovar_path}/convert2annovar.pl \
#	-format vcf4 ${resultdir}/${sample}.filtered.merged.split.DP30.vcf \
#	-outfile ${resultdir}/${sample}.DP30.processed_gnomAD.avinput
#	
#
#python ${gnomad2_1_annotate_oea} \
#	${resultdir}/${sample}.DP30.processed_gnomAD.avinput \
#	${resultdir}/${sample}.DP30.processed_gnomAD2.1.hg19.txt
#
#python ${gnomad2_1_AF} \
#	${resultdir}/${sample}.DP30.processed_gnomAD2.1.hg19.txt \
#	${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.avinput
#
#cp ${intervar_path}/config.ini ./
##
#${intervar_path}/./Intervar.py \
#	-b hg19 \
#	-i ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.avinput \
#	--input_type=AVinput \
#	-t ${intervar_path}/intervardb/ \
#	-d ${intervar_path}/humandb/ \
#	-o ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001
###
##
#python ${filter_intervar_LOF} \
#	${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.hg19_multianno.txt.intervar \
#	${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.LOF.intervar
#
#
#egrep -v -i "InterVar: Benign|InterVar: Likely benign" ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.LOF.intervar > ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.debenign.intervar
#
#egrep -i "InterVar: Pathogenic|InterVar: Likely pathogenic" ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.LOF.intervar > ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.patho.intervar
#
##
#sed 's/^/'"${sample}"'\t&/g' ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.debenign.intervar > ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.debenign.intervar.withid
##
##sed 's/^/'"${sample}"'\t&/g' ${resultdir}/${sample}.sb.af0.2dp10_variants_gnomAD0.001.sb.LOF.intervar > ${resultdir}/${sample}.sb.af0.2dp10_variants_gnomAD0.001.sb.LOF.intervar.withid
#
#sed 's/^/'"${sample}"'\t&/g' ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.patho.intervar > ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.patho.intervar.withid
#
#python ${check_zygo} \
#	-i ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.debenign.intervar.withid \
#	-a ${resultdir}/${sample}.DP30.processed_gnomAD.avinput \
#	-o ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.debenign.intervar.withid.zygosity
#######################################################################
#python ${check_zygo} \
#	-i ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.patho.intervar.withid \
#	-a ${resultdir}/${sample}.DP30.processed_gnomAD.avinput \
#	-o ${resultdir}/${sample}.DP30.processed_gnomAD2.1.0.001.patho.intervar.withid.zygosity
#############HEREIN process GnomAD2.0 results########################

#python ${gnomad2_0_AF}	\
#	${resultdir}/${sample}.processed_gnomAD.hg19_multianno.txt \
#	${resultdir}/${sample}.processed_gnomAD0.001.avinput
#
#cp ${intervar_path}/config.ini ./
##
#${intervar_path}/./Intervar.py \
#	-b hg19 \
#	-i ${resultdir}/${sample}.processed_gnomAD0.001.avinput \
#	--input_type=AVinput \
#	-t ${intervar_path}/intervardb/ \
#	-d ${intervar_path}/humandb/ \
#	-o ${resultdir}/${sample}.processed_gnomAD0.001
##
##
#python ${filter_intervar_LOF} \
#	${resultdir}/${sample}.processed_gnomAD0.001.hg19_multianno.txt.intervar \
#	${resultdir}/${sample}.processed_gnomAD0.001.LOF.intervar
#
#
#egrep -v -i "InterVar: Benign|InterVar: Likely benign" ${resultdir}/${sample}.processed_gnomAD0.001.LOF.intervar > ${resultdir}/${sample}.processed_gnomAD0.001.debenign.intervar
#
#egrep -i "InterVar: Pathogenic|InterVar: Likely pathogenic" ${resultdir}/${sample}.processed_gnomAD0.001.LOF.intervar > ${resultdir}/${sample}.processed_gnomAD0.001.patho.intervar
#
##
#sed 's/^/'"${sample}"'\t&/g' ${resultdir}/${sample}.processed_gnomAD0.001.debenign.intervar > ${resultdir}/${sample}.processed_gnomAD0.001.debenign.intervar.withid
##
##sed 's/^/'"${sample}"'\t&/g' ${resultdir}/${sample}.sb.af0.2dp10_variants_gnomAD0.001.sb.LOF.intervar > ${resultdir}/${sample}.sb.af0.2dp10_variants_gnomAD0.001.sb.LOF.intervar.withid
#
#sed 's/^/'"${sample}"'\t&/g' ${resultdir}/${sample}.processed_gnomAD0.001.patho.intervar > ${resultdir}/${sample}.processed_gnomAD0.001.patho.intervar.withid
##
###python ${ad_dp_filter_dp20} ${resultdir}/${sample}.filtered_indels.vcf ${resultdir}/${sample}.filtered_addp20_indels.vcf
###########InterVar###########
#cp ${intervar_path}/config.ini ./
#
#${intervar_path}/./Intervar.py \
#	-b hg19 \
#	-i ${resultdir}/${sample}.filtered_addp30_snps.vcf \
#	--input_type=VCF \
#	-t ${intervar_path}/intervardb/ \
#	-d ${intervar_path}/humandb/ \
#	-o ${resultdir}/${sample}.snp30.intervar
# 
#${intervar_path}/./Intervar.py \
#	-b hg19 \
#	-i ${resultdir}/${sample}.filtered_addp30_indels.vcf \
#	--input_type=VCF \
#	-t ${intervar_path}/intervardb/ \
#	-d ${intervar_path}/humandb/ \
#	-o ${resultdir}/${sample}.indel30.intervar
#
###########Select-pathogenic-variants################
#cut -f 1-12,14 ${resultdir}/${sample}.snp30.intervar.*.intervar|egrep -i "#chr|pathogenic" > ${resultdir}/${sample}.snp30.patho.intervar
# 
#cut -f 1-12,14 ${resultdir}/${sample}.indel30.intervar.*.intervar|egrep -i "#chr|pathogenic" > ${resultdir}/${sample}.indel30.patho.intervar
#
#cat ${resultdir}/${sample}.snp30.patho.intervar ${resultdir}/${sample}.indel30.patho.intervar > ${resultdir}/${sample}.all30.patho.intervar
#
#grep -v -i "#chr" ${resultdir}/${sample}.snp30.patho.intervar | cut -f 1-2 > ${resultdir}/${sample}.snp30.patho.pos
#
######NOTE:: NEED TO USE THE DEHCTAG FILES TO MAKE THE DP CONSISTENT###
#${samtools} view \
#  -b -r "${sample}" \
#  ${bamdir}/${sample}.sb.final.bam > ${resultdir}/${sample}.deHCtag.bam
#
#${samtools} index \
#  -b -@ "${threads}" \
#  ${resultdir}/${sample}.deHCtag.bam ${resultdir}/${sample}.deHCtag.bai
#
${samtools} mpileup \
  --reference ${ref37} \
  -l ${maindir}/bin/BJCH_debenign_variant.noSB.sort.bed \
  -q 0 -Q 0 -s -O -a -a \
  --output-QNAME --output-MQ \
  -o ${resultdir}/${sample}.noSB.mpileup \
  ${resultdir}/${sample}.deHCtag.bam
#
#${samtools} depth \
#	--reference ${ref37} \
#	-a -a -b ${maindir}/bin/BJCH_debenign_variant.sort.bed \
#	${resultdir}/${sample}.deHCtag.bam > ${resultdir}/${sample}.final.depth

sed 's/^/'"${sample}"'\t&/g' ${resultdir}/${sample}.noSB.mpileup > ${resultdir}/${sample}.noSB.mpileup.withid

#python ${pileup_sum} ${resultdir}/${sample}.snp30.mpileup ${resultdir}/${sample}.snp30.mpileup.sum.txt ${sample}

#${intervar_path}/./Intervar.py \
#	-b hg19 \
#	-i ${resultdir}/${sample}.filtered_addp20_snps.vcf \
#	--input_type=VCF \
#	-t ${intervar_path}/intervardb/ \
#	-d ${intervar_path}/humandb/ \
#	-o ${resultdir}/${sample}.snp20.intervar
# 
#${intervar_path}/./Intervar.py \
#	-b hg19 \
#	-i ${resultdir}/${sample}.filtered_addp20_indels.vcf \
#	--input_type=VCF \
#	-t ${intervar_path}/intervardb/ \
#	-d ${intervar_path}/humandb/ \
#	-o ${resultdir}/${sample}.indel20.intervar
#
###########Select-pathogenic-variants################
#cut -f 1-12,14 ${resultdir}/${sample}.snp20.intervar.*.intervar|egrep -i "#chr|pathogenic" > ${resultdir}/${sample}.snp20.patho.intervar
# 
#cut -f 1-12,14 ${resultdir}/${sample}.indel20.intervar.*.intervar|egrep -i "#chr|pathogenic" > ${resultdir}/${sample}.indel20.patho.intervar
#
#cat ${resultdir}/${sample}.snp20.patho.intervar ${resultdir}/${sample}.indel20.patho.intervar > ${resultdir}/${sample}.all20.patho.intervar
#
#grep -v -i "#chr" ${resultdir}/${sample}.snp20.patho.intervar | cut -f 1-2 > ${resultdir}/${sample}.snp20.patho.pos
#
######NOTE:: NEED TO USE THE DEHCTAG FILES TO MAKE THE DP CONSISTENT###
#
#${samtools} mpileup \
#  --reference ${ref37} \
#  -l ${resultdir}/${sample}.snp20.patho.pos \
#  -q 0 -Q 0 -s -O \
#  --output-QNAME --output-MQ \
#  -o ${resultdir}/${sample}.snp20.mpileup \
#  ${resultdir}/${sample}.deHCtag.bam
#
#python ${pileup_sum} ${resultdir}/${sample}.snp20.mpileup ${resultdir}/${sample}.snp20.mpileup.sum.txt ${sample}

#python ${mpile_sum} ${resultdir}/${sample}.snp20.mpileup.sum.txt ${resultdir}/${sample}.snp20.mpileup.ratio

#python ${mpile_sum} ${resultdir}/${sample}.snp30.mpileup.sum.txt ${resultdir}/${sample}.snp30.mpileup.ratio


#scp ${fqpath}/${sample}.R?.clean.fastq.gz bd2:/home/yanmuqing/BJCH/fastq
#
#rm config.ini
#
#rm ${resultdir}/${sample}.aln.bam
#
#rm ${fqpath}/${sample}.R?.clean.fastq.gz 
