#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
# #SBATCH -J selfmap
#SBATCH -o selfmap_%x_%j.out
#SBATCH -e selfmap_%x_%j.err




source ~/.bashrc

fastq=$1   # $r.palindrome_treated.fq.gz
threads=32

fastq_prefix=$(basename "$fastq") 
fastq_prefix=${fastq_prefix%.fastq.gz}
fastq_prefix=${fastq_prefix%.fasta.gz}
fastq_prefix=${fastq_prefix%.fq.gz}
fastq_prefix=${fastq_prefix%.fa.gz}
fastq_prefix=${fastq_prefix%.fastq}
fastq_prefix=${fastq_prefix%.fasta}
fastq_prefix=${fastq_prefix%.fq}
fastq_prefix=${fastq_prefix%.fa} # $r.ccs

asm_prefix=$fastq_prefix.asm


echo 'processing '$r'...'

dm6_ref=/public5/home/sch5655/dm6.fa

conda activate base
echo 'hifiasm assembling ...'
srun -N1 -n 1 -c $threads hifiasm -t $threads -l0 -o $asm_prefix -f0 $fastq
srun -N1 -n 1 -c 1 awk '/^S/{print ">"$2;print $3}' $asm_prefix.bp.p_ctg.gfa > $asm_prefix.bp.p_ctg.fa
echo 'hifiasm done'

echo 'quast large ...'
srun -N1 -n 1 -c $threads quast.py --large -t $threads -r $dm6_ref -o $asm_prefix.p.quastLG $asm_prefix.bp.p_ctg.fa
echo 'quast done'

echo 'merqury ...'
srun -N 1 -n 1 -c $threads meryl count k=19 output $fastq_prefix.meryl $fastq
srun -N 1 -n 1 -c $threads merqury.sh $fastq_prefix.meryl $asm_prefix.bp.p_ctg.fa $asm_prefix.bp.p_ctg
echo 'merqury done'

echo 'compleasm ...'
conda activate compleasm
srun -N 1 -n 1 -c $threads compleasm run -t $threads -l diptera -a $asm_prefix.bp.p_ctg.fa -o $asm_prefix.bp.p_ctg.compleasm -L /public5/home/sch5655/wmx/compleasm_download/
echo 'compleasm done'





