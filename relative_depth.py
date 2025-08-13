import sys
import pandas as pd

def main():
    indepth = sys.argv[1] # *.depth

    binsize = 100000
    out = ".".join(indepth.split("/")[-1].split(".")[0:-1]) + "." + str(binsize) + ".relative.depth"

    result = mean_depth_per_bin(indepth, binsize)
    median_bindepth = median_depth_without_chrXY(result)
    result['relative_depth'] = result['mean_depth'] / (median_bindepth / 2)

    result = result.reset_index(drop=True)
    result['bin_number'] = result.index + 1

    # Export as requested: chrom, bin_start, relative_depth, bin_number
    result[['chrom', 'bin_start', 'relative_depth', 'bin_number']].to_csv(
        out, sep='\t', index=False, header=False, quoting=3)


def median_depth_without_chrXY(result):
    autosomes = result[~result['chrom'].isin(['chrX', 'chrY'])]
    return autosomes['mean_depth'].median()

def mean_depth_per_bin(indepth, binsize):
    df = pd.read_csv(indepth, sep='\t', header=None, names=['chrom', 'pos', 'depth'])
    df['bin'] = ((df['pos'] - 1) // binsize) * binsize + 1
    result = (
        df.groupby(['chrom', 'bin'])['depth']
        .mean()
        .reset_index()
        .rename(columns={'bin': 'bin_start', 'depth': 'mean_depth'})
    )
    return result


if __name__ == "__main__":
    main()

