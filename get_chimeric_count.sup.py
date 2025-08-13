import sys
import pysam

def main():
    inbam = sys.argv[1] # inbam (-F4 bam)
    out_name = sys.argv[2] # H54 H58

    inbam_dic = read_inbam(inbam)
    chimeric_count, total_count = count_chimeric_sup(inbam_dic)
    # out_line = inbam.split("/")[-1] + "\t" + str(chimeric_count) + "\t" + str(total_count) + "\t" + str(chimeric_count/total_count)
    out_line = out_name + "\tsupplementary_mapping_chimera\t" + str(chimeric_count) + "\t" + str(total_count) + "\t" + str(chimeric_count/total_count)
    print(out_line)


def read_inbam(inbam):
    inbam_dic = {}
    with pysam.AlignmentFile(inbam, "rb") as bamfile:
        for read in bamfile:
            if not read.is_unmapped and not read.is_secondary:
                readname = read.query_name
                ref_pos = read.reference_name + "\t" + str(read.reference_start)
                inbam_dic.setdefault(readname, []).append(ref_pos)
    return inbam_dic

def count_chimeric_sup(inbam_dic):
    total_count = 0
    chimeric_count = 0
    for rname in inbam_dic:
        total_count += 1
        if len(inbam_dic[rname]) > 1:
            chimeric_count += 1
    return chimeric_count, total_count


# def count_chimeric(inbam):
#     chimeric_count = 0
#     total_count = 0
#     with pysam.AlignmentFile(inbam, "rb") as bamfile:
#         for read in bamfile:
#             if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
#                 total_count += 1
#                 read_cigartuples = read.cigartuples
#                 if ((read_cigartuples[0][0] == 4 and read_cigartuples[0][1] > 50) or (read_cigartuples[0][0] == 5 and read_cigartuples[0][1] > 50)) or ((read_cigartuples[-1][0] == 4 and read_cigartuples[-1][1] > 50) or (read_cigartuples[-1][0] == 5 and read_cigartuples[-1][1] > 50)): # Soft clip (clipped sequence present in SEQ) / Hard clip (clipped sequence absent from SEQ)
#                     chimeric_count += 1
#     return chimeric_count, total_count


if __name__ == "__main__":
    main()
