import sys
import pysam

def main():
    in_chimeric_bam = sys.argv[1]
    out = in_chimeric_bam.split("/")[-1] + ".inverted_repeat.rname"

    all_position_dic, primary_mapping_position_dic = read_inbam_get_primary_mapping(in_chimeric_bam)

    chimeric_count, total_count = count_chimeric_sup(all_position_dic)
    inverted_repeat_count = count_inverted_repeat(primary_mapping_position_dic, all_position_dic, out)

    print(str(inverted_repeat_count) + "\t" + str(chimeric_count) + "\t" + str(total_count))


def read_inbam_get_primary_mapping(inbam):
    all_position_dic = {}
    primary_mapping_position_dic = {}
    with pysam.AlignmentFile(inbam, "rb") as bamfile:
        for read in bamfile:
            readname = read.query_name
            chr_name = read.reference_name
            start = read.reference_start
            end = read.reference_end
            if read.is_reverse:
                strand = "-"
            else:
                strand = "+"

            ref_pos = chr_name + "\t" + str(start) + "\t" + str(end) + "\t" + strand
            all_position_dic.setdefault(readname, []).append(ref_pos)

            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                primary_mapping_position_dic[readname] = ref_pos

    return all_position_dic, primary_mapping_position_dic


def count_chimeric_sup(all_position_dic):
    total_count = 0
    chimeric_count = 0
    for rname in all_position_dic:
        total_count += 1
        if len(all_position_dic[rname]) > 1:
            chimeric_count += 1
    return chimeric_count, total_count


def if_inverted_repeat(primary_pos, pos):
    primary_chr, primary_start, primary_end, primary_strand = primary_pos.split("\t")
    pos_chr, pos_start, pos_end, pos_strand = pos.split("\t")
    if primary_chr == pos_chr and pos_start >= primary_start and pos_end <= primary_end:
        if primary_strand != pos_strand:
            return True
    return False


def count_inverted_repeat(primary_mapping_position_dic, all_position_dic, out):
    with open(out, "w") as final:
        inverted_repeat_count = 0
        for rname in primary_mapping_position_dic:
            primary_pos = primary_mapping_position_dic[rname]
            if primary_pos in all_position_dic[rname] and len(all_position_dic[rname]) > 1:
                for pos in all_position_dic[rname]:
                    if pos != primary_pos:
                        if if_inverted_repeat(primary_pos, pos):
                            inverted_repeat_count += 1
                            final.write(rname + "\n")
                            break
    return inverted_repeat_count


if __name__ == "__main__":
    main()
