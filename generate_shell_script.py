import sys
import os

def main():
    intsv = sys.argv[1]
    jobname = sys.argv[2]

    total_cores = 64

    out = "a.sh"

    original_ccs_bam_prefix_dic, sample_count, sample_ls = read_tsv(intsv)
    generated_out_shell_script(original_ccs_bam_prefix_dic, total_cores, sample_count, sample_ls, jobname, out)


def read_tsv(intsv):
    original_ccs_bam_prefix_dic = {}
    sample_count = 0
    sample_ls = []
    with open(intsv, "r") as inf:
        for line in inf:
            sample_count += 1
            line = line.rstrip("\n")
            ll = line.split("\t")
            prefix = ll[1]

            sample_ls.append(prefix)

            original_ccs_bamfile_absolute_path = ll[0]
            original_ccs_bam_prefix_dic[prefix] = original_ccs_bamfile_absolute_path

    sample_ls.sort()
    return original_ccs_bam_prefix_dic, sample_count, sample_ls

def generated_out_shell_script(original_ccs_bam_prefix_dic, total_cores, sample_count, sample_ls, jobname, out):
    with open(out, "w") as outf:
        header = "#!/bin/bash\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c " + str(total_cores) + "\n#SBATCH -J " + jobname + "\n#SBATCH -o job.out\n#SBATCH -e job.err\n"
        work_dir = "cd " + os.getcwd() + "\n"

        bam2fastq_parts = bam2fastq(original_ccs_bam_prefix_dic, total_cores, sample_count)
        primary_asm_parts = genome_asm(sample_ls)
        self_mapping_palindrom_parts = self_mapping_palindrom_asm(sample_ls)
        mapping_parts = reads_mapping_to_dm6(sample_ls, total_cores)
        chimeric_rate_inverted_repeat_parts = chimeric_rate_inverted_repeat(jobname, sample_ls, total_cores, sample_count)
        relative_depth_parts = relative_depth(sample_ls, total_cores, sample_count)
        data_size_readlength_parts = data_size(jobname, total_cores, sample_count, sample_ls)
        mapping_rate_parts = mapping_rate(jobname, total_cores, sample_count, sample_ls)

        outf.write(header + "\n\n")
        outf.write(work_dir + "\n\n")
        outf.write(bam2fastq_parts + "\n\n")
        outf.write(primary_asm_parts + "\n\n")
        outf.write(self_mapping_palindrom_parts + "\n\n")
        outf.write(mapping_parts + "\n\n")
        outf.write(chimeric_rate_inverted_repeat_parts + "\n\n")
        outf.write(relative_depth_parts + "\n\n")
        outf.write(data_size_readlength_parts + "\n\n")
        outf.write(mapping_rate_parts + "\n\n")



def bam2fastq(original_ccs_bam_prefix_dic, total_cores, sample_count):
    outline = '#### bam2fastq\necho "bam2fastq ..."\n'
    
    per_sample_cores = total_cores // sample_count
    for prefix in dict(sorted(original_ccs_bam_prefix_dic.items())):
        outline0 = "srun -N 1 -n 1 -c " + str(per_sample_cores) + " bam2fastq " + original_ccs_bam_prefix_dic[prefix] + " -o " + prefix + ".ccs &\n"
        outline = outline + outline0
    
    outline = outline + "wait\n" + 'echo "bam2fastq done"\n'
    return outline

def genome_asm(sample_ls):
    outline = '#### primary asm\necho "submitting asm.sh..."\n'

    outline = outline + "for r in " + " ".join(sample_ls) + "\ndo\n"
    outline = outline + '\tsbatch genome_asm_pipeline.sh $r.ccs.fastq.gz\ndone\n'
    outline = outline + 'echo "submitting asm.sh done"\n'
    return outline

def self_mapping_palindrom_asm(sample_ls):
    outline = '#### self mapping + palindrom_find asm\n'
    outline = outline + "for r in " + " ".join(sample_ls) + "\ndo\n"
    outline = outline + '\techo "submitting $r slef_asm.sh..."\n'
    outline = outline + '\tsbatch self_mapping_palindrom_asm.sh $r.ccs.fastq.gz\n'
    outline = outline + '\techo "submitting $r slef_asm.sh done"\ndone\n'
    return outline 

def reads_mapping_to_dm6(sample_ls, total_cores):
    outline = '#### reads mapping to dm6\necho "Starting mapping to dm6..."\n'
    outline = outline + "for r in " + " ".join(sample_ls) + "\ndo\n"
    outline = outline + '\tdm6_ref=/public5/home/sch5655/dm6.fa'
    outline = outline + '\tfastq=$r.ccs.fastq.gz\n\n'
    outline = outline + '\t### mapping\n'
    outline = outline + '\tsrun -N 1 -n 1 -c 64 minimap2 -t ' + str(total_cores) + ' -ax map-hifi --MD --secondary=no $dm6_ref $fastq | \\\n'
    outline = outline + '\tsamtools view -@' + str(total_cores) + ' -bS -F0x4 > $r.ccs.fq.dm6.F4.bam\n'
    outline = outline + '\tsrun -N 1 -n 1 -c ' + str(total_cores) + ' samtools sort -@' + str(total_cores) + ' $r.ccs.fq.dm6.F4.bam -o $r.ccs.fq.dm6.F4.s.bam\n'
    outline = outline + '\tsrun -N 1 -n 1 -c ' + str(total_cores) + ' samtools index -@' + str(total_cores) + ' $r.ccs.fq.dm6.F4.s.bam\n'
    outline = outline + '\trm $r.ccs.fq.dm6.F4.bam\n\n'
    outline = outline + '\tsamtools depth -@' + str(total_cores) + ' -a $r.ccs.fq.dm6.F4.s.bam | awk \'{if($1 !~ /_/ && $1 !~ /chrM/){print $0}}\' > $r.ccs.fq.dm6.F4.s.bam.depth\n\n'
    outline = outline + 'done\n'
    outline = outline + 'echo "Starting mapping to dm6 done"\n'
    return outline

def chimeric_rate_inverted_repeat(jobname, sample_ls, total_cores, sample_count):
    outline = '#### chimeric rate and inverted repeat count\necho "Starting calculating chimeric rate and inverted repeat count..."\n'
    outline = outline + '> ' + str(jobname) + '.ccs.summary.chimericrate.sup.tsv\n'
    outline = outline + '> ' + str(jobname) + '.ccs.inverted_repeat.count\n'

    per_job_cores = total_cores // (sample_count * 2)

    for r in sample_ls:
        chimeric_rate_line = "srun -N 1 -n 1 -c " + str(per_job_cores) + " python3 /public5/home/sch5655/wmx/get_chimeric_count.sup.py " + r + ".ccs.fq.dm6.F4.s.bam " + r + " >> " + str(jobname) + ".ccs.summary.chimericrate.sup.tsv &\n"
        inverted_reepat_line = "srun -N 1 -n 1 -c " + str(per_job_cores) + " python3 /public5/home/sch5655/wmx/get_inverted_repeat.py " + r + ".ccs.fq.dm6.F4.s.bam | awk '{print \"" + r + "\\t\"\$0}' >> " + str(jobname) + ".ccs.inverted_repeat.count &\n"
        outline = outline + chimeric_rate_line + inverted_reepat_line
    outline = outline + 'wait\necho "Starting calculating chimeric rate and inverted repeat count done"\n'

    return outline

def relative_depth(sample_ls, total_cores, sample_count):
    outline = '#### relative depth\necho "Starting calculating relative depth..."\n'
    outline = outline + '#### relative depth\n'

    per_job_cores = total_cores // sample_count
    for r in sample_ls:
        relative_depth_line = "srun -N 1 -n 1 -c " + str(per_job_cores) + " python3 /public5/home/sch5655/wmx/relative_depth.py " + r + ".ccs.fq.dm6.F4.s.bam.depth &\n"
        outline = outline + relative_depth_line
    outline = outline + 'wait\nnecho "Starting calculating relative depth done"\n'

    return outline


def data_size(jobname, total_cores, sample_count, sample_ls):
    outline = '#### data size\n'
    outline = outline + '>' + str(jobname) + '.samples.total_basecount.out\n'    

    per_job_cores = total_cores // sample_count
    for r in sample_ls:
        relative_depth_line = "srun -N 1 -n 1 -c " + str(per_job_cores) + " pigz -dc " + r + ".ccs.fastq.gz | awk 'BEGIN{total_len = 0}{if(NR % 4 == 2){total_len = total_len + length($0)}}END{print \"" + r + "\\t\"total_len}' >> " + str(jobname) + ".samples.total_basecount.out &\n"
        outline = outline + relative_depth_line
    outline = outline + 'wait\n\n'

    for r in sample_ls:
        readlength_line = "srun -N 1 -n 1 -c " + str(per_job_cores) + " pigz -dc " + r + ".ccs.fastq.gz | awk '{if(NR % 4 == 2){print length($0)}}' > " + r + ".ccs.fastq.readlength &\n"
        outline = outline + readlength_line
    outline = outline + 'wait\n'

    return outline


def mapping_rate(jobname, total_cores, sample_count, sample_ls):
    outline = '#### mapping rate\n'
    outline = outline + '> ' + str(jobname) + '.ccs.fastq.readcount\n'

    per_job_cores = total_cores // sample_count
    for r in sample_ls:
        readcount_line = "srun -N 1 -n 1 -c " + str(per_job_cores) + "awk 'END{print \"" + r + "\\t\"NR}' " + r + ".ccs.fastq.readlength >> " + str(jobname) + ".ccs.fastq.readcount &\n"
        outline = outline + readcount_line
    outline = outline + 'wait\n\n'

    outline = outline + '> ' + str(jobname) + '.ccs.fastq.mapping.readcount\n'
    outline = outline + "for r in " + " ".join(sample_ls) + "\ndo\n"
    outline = outline + '\tinbam=$r.ccs.fq.dm6.F4.s.bam\n\n'
    outline = outline + '\tsamtools view -@' + str(total_cores) + ' -F0x904 $inbam | awk \'BEGIN{c = 0}{c++}END{print "$r\\t"c}\' >> ' + str(jobname) + '.ccs.fastq.mapping.readcount &\n'
    outline = outline + 'done\n\n'
    outline = outline + "awk '{if(NR == FNR){all[$1] = $2} else{if($1 in all){rate = $2 / all[$1]; print $1\"\\t\"$2\"\\t\"all[$1]\"\\t\"rate}}}' " + str(jobname) + ".ccs.fastq.readcount " + str(jobname) + ".ccs.fastq.mapping.readcount > " + str(jobname) + ".ccs.fastq.mapping.rate\n"

    return outline


if __name__ == "__main__":
    main()


