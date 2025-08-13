import sys
import pandas as pd

def main():
    indepth = sys.argv[1] # *.depth

    binsize = 100000
    out = ".".join(indepth.split("/")[-1].split(".")[0:-1]) + "." + str(binsize) + ".relative.depth1"

    bin_mean_depth_dic = mean_depth_per_bin_stream(indepth, binsize)
    export_bin_depth1(bin_mean_depth_dic, out)


def mean_depth_per_bin_stream(indepth, binsize):
    bin_id = 0
    bin_mean_depth_dic = {}
    with open(indepth, 'r') as f:
        for line in f:
            ll = line.rstrip("\n").split("\t")
            chrom = ll[0]
            pos = int(ll[1])
            depth = int(ll[2])
            bin_start = ((pos - 1) // binsize) * binsize + 1

            if chrom in bin_mean_depth_dic and bin_start in bin_mean_depth_dic[chrom]:
                bin_mean_depth_dic[chrom][bin_start][0] += depth
                bin_mean_depth_dic[chrom][bin_start][1] += 1
            else:
                bin_id += 1
                bin_mean_depth_dic.setdefault(chrom, {})[bin_start] = [depth, 1, bin_id]
    return bin_mean_depth_dic

def export_bin_depth1(bin_mean_depth_dic, out):
    with open(out, "w") as outf:
        for chrom in bin_mean_depth_dic:
            for bin_start in bin_mean_depth_dic[chrom]:
                depth, count, bin_id = bin_mean_depth_dic[chrom][bin_start]
                mean_depth = depth / count
                outf.write(f"{chrom}\t{bin_start}\t{mean_depth}\t{bin_id}\n")


if __name__ == "__main__":
    main()
