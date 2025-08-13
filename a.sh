#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -J pta_0724
#SBATCH -o pta.out
#SBATCH -e pta.err



cd /public5/home/sch5655/wmx/PTA_20250724/data_process


### 2025.7.17
echo "bam2fastq ..."
srun -N 1 -n 1 -c 8 bam2fastq /public5/home/sch5655/wmx/PTA_20250724/25TW01557/r84258_250716_001_1_C01/XT1-1.bc2065.bc2065.HiFi.bam -o XT1_1.ccs &
srun -N 1 -n 1 -c 8 bam2fastq /public5/home/sch5655/wmx/PTA_20250724/25TW01558/r84258_250716_001_1_C01/XT1-0.bc2066.bc2066.HiFi.bam -o XT1_0.ccs &
srun -N 1 -n 1 -c 8 bam2fastq /public5/home/sch5655/wmx/PTA_20250724/25TW01559/r84258_250716_001_1_C01/XT2-1.bc2067.bc2067.HiFi.bam -o XT2_1.ccs &
srun -N 1 -n 1 -c 8 bam2fastq /public5/home/sch5655/wmx/PTA_20250724/25TW01560/r84258_250716_001_1_C01/XT2-0.bc2068.bc2068.HiFi.bam -o XT2_0.ccs &
srun -N 1 -n 1 -c 8 bam2fastq /public5/home/sch5655/wmx/PTA_20250724/25TW01561/r84130_250712_001_1_C01/E1-1.bc2069.bc2069.HiFi.bam -o E1_1.ccs &
srun -N 1 -n 1 -c 8 bam2fastq /public5/home/sch5655/wmx/PTA_20250724/25TW01562/r84130_250712_001_1_C01/E1-0.bc2070.bc2070.HiFi.bam -o E1_0.ccs &
srun -N 1 -n 1 -c 8 bam2fastq /public5/home/sch5655/wmx/PTA_20250724/25TW01563/r84130_250712_001_1_C01/E2-1.bc2071.bc2071.HiFi.bam -o E2_1.ccs &
srun -N 1 -n 1 -c 8 bam2fastq /public5/home/sch5655/wmx/PTA_20250724/25TW01564/r84130_250712_001_1_C01/E2-0.bc2072.bc2072.HiFi.bam -o E2_0.ccs &
wait
echo "bam2fastq done"



#### primary asm
echo "submitting asm.sh..."
sbatch asm.sh
echo "submitting asm.sh done"


#### self mapping + palindrom_find asm
for r in XT1_1 XT1_0 XT2_1 XT2_0 E1_1 E1_0 E2_1 E2_0
do
	echo "submitting $r slef_asm.sh..."
    
	sbatch -J $r.palindrom slef_asm.sh $r

	echo "submitting $r asm.sh done"
done



echo "Starting mapping to dm6..."
for r in XT1_1 XT1_0 XT2_1 XT2_0 E1_1 E1_0 E2_1 E2_0
do
	dm6_ref=/public5/home/sch5655/dm6.fa
	fastq=$r.ccs.fastq.gz

	### mapping
	srun -N 1 -n 1 -c 64 minimap2 -t 64 -ax map-hifi --MD --secondary=no $dm6_ref $fastq | \
	samtools view -@64 -bS -F0x4 > $r.ccs.fq.dm6.F4.bam
	srun -N 1 -n 1 -c 64 samtools sort -@64 $r.ccs.fq.dm6.F4.bam -o $r.ccs.fq.dm6.F4.s.bam
	srun -N 1 -n 1 -c 64 samtools index -@64 $r.ccs.fq.dm6.F4.s.bam
	rm $r.ccs.fq.dm6.F4.bam

	samtools depth -@64 -a $r.ccs.fq.dm6.F4.s.bam | awk '{if($1 !~ /_/ && $1 !~ /chrM/){print $0}}' > $r.ccs.fq.dm6.F4.s.bam.depth

done
echo "Starting mapping to dm6 done"





echo "Starting calculating chimeric rate and inverted repeat count..."
> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.summary.chimericrate.sup.tsv
> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.inverted_repeat.count
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_chimeric_count.sup.py XT1_1.ccs.fq.dm6.F4.s.bam XT1_1 >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.summary.chimericrate.sup.tsv &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_chimeric_count.sup.py XT1_0.ccs.fq.dm6.F4.s.bam XT1_0 >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.summary.chimericrate.sup.tsv &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_chimeric_count.sup.py XT2_1.ccs.fq.dm6.F4.s.bam XT2_1 >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.summary.chimericrate.sup.tsv &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_chimeric_count.sup.py XT2_0.ccs.fq.dm6.F4.s.bam XT2_0 >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.summary.chimericrate.sup.tsv &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_chimeric_count.sup.py E1_1.ccs.fq.dm6.F4.s.bam E1_1 >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.summary.chimericrate.sup.tsv &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_chimeric_count.sup.py E1_0.ccs.fq.dm6.F4.s.bam E1_0 >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.summary.chimericrate.sup.tsv &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_chimeric_count.sup.py E2_1.ccs.fq.dm6.F4.s.bam E2_1 >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.summary.chimericrate.sup.tsv &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_chimeric_count.sup.py E2_0.ccs.fq.dm6.F4.s.bam E2_0 >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.summary.chimericrate.sup.tsv &

srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_inverted_repeat.py XT1_1.ccs.fq.dm6.F4.s.bam | awk '{print "XT1_1\t"$0}' >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.inverted_repeat.count &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_inverted_repeat.py XT1_0.ccs.fq.dm6.F4.s.bam | awk '{print "XT1_0\t"$0}' >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.inverted_repeat.count &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_inverted_repeat.py XT2_1.ccs.fq.dm6.F4.s.bam | awk '{print "XT2_1\t"$0}' >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.inverted_repeat.count &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_inverted_repeat.py XT2_0.ccs.fq.dm6.F4.s.bam | awk '{print "XT2_0\t"$0}' >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.inverted_repeat.count &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_inverted_repeat.py E1_1.ccs.fq.dm6.F4.s.bam | awk '{print "E1_1\t"$0}' >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.inverted_repeat.count &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_inverted_repeat.py E1_0.ccs.fq.dm6.F4.s.bam | awk '{print "E1_0\t"$0}' >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.inverted_repeat.count &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_inverted_repeat.py E2_1.ccs.fq.dm6.F4.s.bam | awk '{print "E2_1\t"$0}' >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.inverted_repeat.count &
srun -N 1 -n 1 -c 4 python3 /public5/home/sch5655/wmx/get_inverted_repeat.py E2_0.ccs.fq.dm6.F4.s.bam | awk '{print "E2_0\t"$0}' >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.inverted_repeat.count &
wait


#### relative depth
srun -N 1 -n 1 -c 8 python3 /public5/home/sch5655/wmx/relative_depth.py XT1_1.ccs.fq.dm6.F4.s.bam.depth &
srun -N 1 -n 1 -c 8 python3 /public5/home/sch5655/wmx/relative_depth.py XT1_0.ccs.fq.dm6.F4.s.bam.depth &
srun -N 1 -n 1 -c 8 python3 /public5/home/sch5655/wmx/relative_depth.py XT2_1.ccs.fq.dm6.F4.s.bam.depth &
srun -N 1 -n 1 -c 8 python3 /public5/home/sch5655/wmx/relative_depth.py XT2_0.ccs.fq.dm6.F4.s.bam.depth &
srun -N 1 -n 1 -c 8 python3 /public5/home/sch5655/wmx/relative_depth.py E1_1.ccs.fq.dm6.F4.s.bam.depth &
srun -N 1 -n 1 -c 8 python3 /public5/home/sch5655/wmx/relative_depth.py E1_0.ccs.fq.dm6.F4.s.bam.depth &
srun -N 1 -n 1 -c 8 python3 /public5/home/sch5655/wmx/relative_depth.py E2_1.ccs.fq.dm6.F4.s.bam.depth &
srun -N 1 -n 1 -c 8 python3 /public5/home/sch5655/wmx/relative_depth.py E2_0.ccs.fq.dm6.F4.s.bam.depth &
wait

echo "Starting calculating chimeric rate and inverted repeat count done"




### 2025.7.25
#### data size
>samples.total_basecount.out
srun -N 1 -n 1 -c 8 pigz -dc E1_0.ccs.fastq.gz | awk 'BEGIN{total_len = 0}{if(NR % 4 == 2){total_len = total_len + length($0)}}END{print "E1_0\t"total_len}' >> samples.total_basecount.out &
srun -N 1 -n 1 -c 8 pigz -dc E1_1.ccs.fastq.gz | awk 'BEGIN{total_len = 0}{if(NR % 4 == 2){total_len = total_len + length($0)}}END{print "E1_1\t"total_len}' >> samples.total_basecount.out &
srun -N 1 -n 1 -c 8 pigz -dc E2_0.ccs.fastq.gz | awk 'BEGIN{total_len = 0}{if(NR % 4 == 2){total_len = total_len + length($0)}}END{print "E2_0\t"total_len}' >> samples.total_basecount.out &
srun -N 1 -n 1 -c 8 pigz -dc E2_1.ccs.fastq.gz | awk 'BEGIN{total_len = 0}{if(NR % 4 == 2){total_len = total_len + length($0)}}END{print "E2_1\t"total_len}' >> samples.total_basecount.out &
srun -N 1 -n 1 -c 8 pigz -dc XT1_0.ccs.fastq.gz | awk 'BEGIN{total_len = 0}{if(NR % 4 == 2){total_len = total_len + length($0)}}END{print "XT1_0\t"total_len}' >> samples.total_basecount.out &
srun -N 1 -n 1 -c 8 pigz -dc XT1_1.ccs.fastq.gz | awk 'BEGIN{total_len = 0}{if(NR % 4 == 2){total_len = total_len + length($0)}}END{print "XT1_1\t"total_len}' >> samples.total_basecount.out &
srun -N 1 -n 1 -c 8 pigz -dc XT2_0.ccs.fastq.gz | awk 'BEGIN{total_len = 0}{if(NR % 4 == 2){total_len = total_len + length($0)}}END{print "XT2_0\t"total_len}' >> samples.total_basecount.out &
srun -N 1 -n 1 -c 8 pigz -dc XT2_1.ccs.fastq.gz | awk 'BEGIN{total_len = 0}{if(NR % 4 == 2){total_len = total_len + length($0)}}END{print "XT2_1\t"total_len}' >> samples.total_basecount.out &
wait

srun -N 1 -n 1 -c 8 pigz -dc E1_0.ccs.fastq.gz | awk '{if(NR % 4 == 2){print length($0)}}' > E1_0.ccs.fastq.readlength &
srun -N 1 -n 1 -c 8 pigz -dc E1_1.ccs.fastq.gz | awk '{if(NR % 4 == 2){print length($0)}}' > E1_1.ccs.fastq.readlength &
srun -N 1 -n 1 -c 8 pigz -dc E2_0.ccs.fastq.gz | awk '{if(NR % 4 == 2){print length($0)}}' > E2_0.ccs.fastq.readlength &
srun -N 1 -n 1 -c 8 pigz -dc E2_1.ccs.fastq.gz | awk '{if(NR % 4 == 2){print length($0)}}' > E2_1.ccs.fastq.readlength &
srun -N 1 -n 1 -c 8 pigz -dc XT1_0.ccs.fastq.gz | awk '{if(NR % 4 == 2){print length($0)}}' > XT1_0.ccs.fastq.readlength &
srun -N 1 -n 1 -c 8 pigz -dc XT1_1.ccs.fastq.gz | awk '{if(NR % 4 == 2){print length($0)}}' > XT1_1.ccs.fastq.readlength &
srun -N 1 -n 1 -c 8 pigz -dc XT2_0.ccs.fastq.gz | awk '{if(NR % 4 == 2){print length($0)}}' > XT2_0.ccs.fastq.readlength &
srun -N 1 -n 1 -c 8 pigz -dc XT2_1.ccs.fastq.gz | awk '{if(NR % 4 == 2){print length($0)}}' > XT2_1.ccs.fastq.readlength &
wait




### 2025.7.31
source ~/.bashrc

dm6_ref=/public5/home/sch5655/dm6.fa
fastq=E2_1.palindrome_treated.fq.gz
asm_prefix=E2_1.palindrome_treated.asm

echo 'processing E2_1...'

conda activate base
echo 'hifiasm assembling ...'
srun -N1 -n 1 -c 64 hifiasm -t 64 -l0 -o $asm_prefix -f0 $fastq
srun -N1 -n 1 -c 1 awk '/^S/{print ">"$2;print $3}' $asm_prefix.bp.p_ctg.gfa > $asm_prefix.bp.p_ctg.fa
echo 'hifiasm done'

echo 'quast large ...'
srun -N1 -n 1 -c 64 quast.py --large -t 64 -r $dm6_ref -o $asm_prefix.p.quastLG $asm_prefix.bp.p_ctg.fa
echo 'quast done'

echo 'merqury ...'
srun -N 1 -n 1 -c 64 meryl count k=19 output E2_1.palindrome_treated.meryl $fastq
srun -N 1 -n 1 -c 64 merqury.sh E2_1.palindrome_treated.meryl $asm_prefix.bp.p_ctg.fa $asm_prefix.bp.p_ctg
echo 'merqury done'

echo 'compleasm ...'
conda activate compleasm
srun -N 1 -n 1 -c 64 compleasm run -t 64 -l diptera -a $asm_prefix.bp.p_ctg.fa -o $asm_prefix.bp.p_ctg.compleasm -L /public5/home/sch5655/wmx/compleasm_download/
echo 'compleasm done'




### 2025.8.13
#### mapping rate
> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.readcount
srun -N 1 -n 1 -c 8 awk 'END{print "E1_0\t"NR}' E1_0.ccs.fastq.readlength >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.readcount &
srun -N 1 -n 1 -c 8 awk 'END{print "E1_0\t"NR}' E1_1.ccs.fastq.readlength >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.readcount &
srun -N 1 -n 1 -c 8 awk 'END{print "E1_0\t"NR}' E2_0.ccs.fastq.readlength >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.readcount &
srun -N 1 -n 1 -c 8 awk 'END{print "E1_0\t"NR}' E2_1.ccs.fastq.readlength >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.readcount &
srun -N 1 -n 1 -c 8 awk 'END{print "E1_0\t"NR}' XT1_0.ccs.fastq.readlength >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.readcount &
srun -N 1 -n 1 -c 8 awk 'END{print "E1_0\t"NR}' XT1_1.ccs.fastq.readlength >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.readcount &
srun -N 1 -n 1 -c 8 awk 'END{print "E1_0\t"NR}' XT2_0.ccs.fastq.readlength >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.readcount &
srun -N 1 -n 1 -c 8 awk 'END{print "E1_0\t"NR}' XT2_1.ccs.fastq.readlength >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.readcount &
wait


>XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.mapping.readcount
for i in XT1_1 XT1_0 XT2_1 XT2_0
do
	inbam=$r.ccs.fq.dm6.F4.s.bam

	samtools view -@64 -F0x904 $inbam | awk 'BEGIN{c = 0}{if($1 in all){a = 1} else{c++; all[$1] = 1}}END{print "'$r'""\t"c}' >> XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.mapping.readcount

done

awk '{if(NR == FNR){all[$1] = $2} else{if($1 in all){rate = $2 / all[$i]; print $1"\t"$2"\t"all[$i]"\t"rate}}}' XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.readcount XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.mapping.readcount > XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.fastq.mapping.rate


