import sys

def main():
    indepth1 = sys.argv[1] # *.depth1

    out_depth2 = ".".join(indepth1.split("/")[-1].split(".")[0:-1]) + ".depth2"

    bin_median_half, all_bin_mean_depth_dic = read_bin_depth1(indepth1)
    export_bin_relative_depth(bin_median_half, all_bin_mean_depth_dic, out_depth2)


def read_bin_depth1(indepth1):
    median_depth_ls = []
    all_bin_mean_depth_dic = {}
    with open(indepth1, "r") as inf:
        for line in inf:
            ll = line.rstrip("\n").split("\t")
            chrom = ll[0]
            bin_start = int(ll[1])
            mean_depth = float(ll[2])
            bin_number = int(ll[3])

            all_bin_mean_depth_dic.setdefault(chrom, {})[bin_start] = [mean_depth, bin_number]

            if chrom != "chrX" and chrom != "chrY":
                median_depth_ls.append(mean_depth)
    
    bin_median_half = sorted(median_depth_ls)[len(median_depth_ls) // 2] / 2
    return bin_median_half, all_bin_mean_depth_dic


def export_bin_relative_depth(bin_median_half, all_bin_mean_depth_dic, out):
    with open(out, "w") as outf:
        for chrom in all_bin_mean_depth_dic:
            for bin_start in all_bin_mean_depth_dic[chrom]:
                mean_depth, bin_number = all_bin_mean_depth_dic[chrom][bin_start]
                relative_depth = mean_depth / bin_median_half
                outf.write(f"{chrom}\t{bin_start}\t{relative_depth}\t{str(bin_number)}\n")

if __name__ == "__main__":
    main()
