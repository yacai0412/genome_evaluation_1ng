#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -J pta_0717_asm
#SBATCH -o asm.out
#SBATCH -e asm.err


cd /public5/home/sch5655/wmx/PTA_20250724/data_process


### 2025.6.3
source ~/.bashrc
for r in XT1_1 XT1_0 XT2_1 XT2_0 E1_1 E1_0 E2_1 E2_0
do
	dm6_ref=/public5/home/sch5655/dm6.fa
	fastq=$r.ccs.fastq.gz

    echo 'processing '$r'...'

	conda activate base
    echo 'hifiasm assembling ...'
	srun -N1 -n 1 -c 64 hifiasm -t 64 -l0 -o $r.ccs.asm -f0 $fastq
	srun -N1 -n 1 -c 1 awk '/^S/{print ">"$2;print $3}' $r.ccs.asm.bp.p_ctg.gfa > $r.ccs.asm.bp.p_ctg.fa
    echo 'hifiasm done'

    echo 'quast large ...'
	srun -N1 -n 1 -c 64 quast.py --large -t 64 -r $dm6_ref -o $r.ccs.asm.p.quastLG $r.ccs.asm.bp.p_ctg.fa
    echo 'quast done'

    echo 'merqury ...'
	srun -N 1 -n 1 -c 64 meryl count k=19 output $r.ccs.meryl $fastq
	srun -N 1 -n 1 -c 64 merqury.sh $r.ccs.meryl $r.ccs.asm.bp.p_ctg.fa $r.ccs.asm.bp.p_ctg
    echo 'merqury done'

	echo 'compleasm ...'
	conda activate compleasm
	srun -N 1 -n 1 -c 64 compleasm run -t 64 -l diptera -a $r.ccs.asm.bp.p_ctg.fa -o $r.ccs.asm.bp.p_ctg.compleasm -L /public5/home/sch5655/wmx/compleasm_download/
    echo 'compleasm done'

done




