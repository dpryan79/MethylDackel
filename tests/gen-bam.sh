# this should have no methylation changes
bwa index ct100.fa
bwa mem ct100.fa ct_R1.fq ct_R2.fq | samtools view -bS - > ct_aln.bam
samtools index ct_aln.bam

bwa index cg100.fa
bwa mem cg100.fa cg_R1.fq cg_R2.fq | samtools view -bS - > cg_aln.bam
samtools index cg_aln.bam
