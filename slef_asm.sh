#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
# #SBATCH -J selfmap
#SBATCH -o selfmap_%x_%j.out
#SBATCH -e selfmap_%x_%j.err

ID0=$1
CPU=32
SAMPLE=${ID0}.ccs.fastq.gz

# srun -N 1 -n 1 -c 32 minimap2 -t $CPU -x ava-pb $SAMPLE $SAMPLE | cat | perl pafIdentifyPalimdrom.pl 1> $ID0.palimProp.1stIte.list 2> $ID0.1stIte.log
# srun -N 1 -n 1 -c 32 perl fastq_partition.and.chop.palindrome.pl $ID0.palimProp.1stIte.list $SAMPLE 1000
# srun -N 1 -n 1 -c 32 minimap2 -t $CPU -x ava-pb $SAMPLE.exclude.fq.gz $SAMPLE.exclude.fq.gz | cat | perl pafIdentifyPalimdrom.pl 1> $ID0.palimProp.2ndIte.list 2> $ID0.2ndIte.log
# srun -N 1 -n 1 -c 32 perl fastq_partition.and.chop.palindrome.pl $ID0.palimProp.2ndIte.list $SAMPLE.exclude.fq.gz 1000
# srun -N 1 -n 1 -c 32 cat $SAMPLE.include.fq.gz $SAMPLE.exclude.fq.gz.exclude.fq.gz $SAMPLE.exclude.fq.gz.include.fq.gz > $ID0.palindrome_treated.fq.gz
# srun -N 1 -n 1 -c 32 rm $SAMPLE.include.fq.gz $SAMPLE.exclude.fq.gz.exclude.fq.gz $SAMPLE.exclude.fq.gz.include.fq.gz $SAMPLE.exclude.fq.gz


#### asm for palindrom spliter
source ~/.bashrc

dm6_ref=/public5/home/sch5655/dm6.fa
fastq=$ID0.palindrome_treated.fq.gz
asm_prefix=$ID0.palindrome_treated.asm

echo 'processing '$ID0'...'

conda activate base
echo 'hifiasm assembling ...'
srun -N1 -n 1 -c 32 hifiasm -t 32 -l0 -o $asm_prefix -f0 $fastq
srun -N1 -n 1 -c 1 awk '/^S/{print ">"$2;print $3}' $asm_prefix.bp.p_ctg.gfa > $asm_prefix.bp.p_ctg.fa
echo 'hifiasm done'

echo 'quast large ...'
srun -N1 -n 1 -c 32 quast.py --large -t 32 -r $dm6_ref -o $asm_prefix.p.quastLG $asm_prefix.bp.p_ctg.fa
echo 'quast done'

echo 'merqury ...'
srun -N 1 -n 1 -c 32 meryl count k=19 output $ID0.palindrome_treated.meryl $fastq
srun -N 1 -n 1 -c 32 merqury.sh $ID0.palindrome_treated.meryl $asm_prefix.bp.p_ctg.fa $asm_prefix.bp.p_ctg
echo 'merqury done'

echo 'compleasm ...'
conda activate compleasm
srun -N 1 -n 1 -c 32 compleasm run -t 32 -l diptera -a $asm_prefix.bp.p_ctg.fa -o $asm_prefix.bp.p_ctg.compleasm -L /public5/home/sch5655/wmx/compleasm_download/
echo 'compleasm done'

