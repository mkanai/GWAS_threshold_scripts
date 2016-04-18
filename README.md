# Empirical estimation of genome-wide significance thresholds based on the 1000 Genomes Project.

This repository contains the scripts used in Kanai, M. *et al.*

## Requirements
* [Python](https://www.python.org/)
* [R](https://www.r-project.org/)
* [GNU Parallel](http://www.gnu.org/software/parallel/)
* Sun Grid Engine (or equivalent) environment
    * We use `qsub` command in our scripts. For other environments, please replace the command appropriately.
* Our customized version of [PLINK 1.9](https://www.cog-genomics.org/plink2/)
    * We added several features to conduct our GWAS simulations against randomly assigned case-control phenotypes.
    * Our customized version is based on [`PLINK v1.90b3l`](https://github.com/chrchang/plink-ng/commit/1d8cbd48106565d381a19efc324472ce47f92e0c).
    * To install, run the following commands.

```{bash}
git clone -b logistic_randomization https://github.com/mkanai/plink-ng
# edit Makefile for your environment
./plink_first_compile
# create a symbolic link
ln -s /path/to/plink-ng/plink /path/to/this/repo/plinklr
```

### R packages
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
* [pforeach](https://github.com/hoxo-m/pforeach)
* [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
* [WRShd](https://github.com/mkanai/WRShd)

To install these packages, run the following commands in R.
```{r}
pkgs <- c("devtools", "data.table", "ggplot2", "Rcpp")
required <- setdiff(pkgs, rownames(installed.packages()))
if (length(required) > 0) {
  install.packages(required)
}
devtools::install_github(c("hoxo-m/pforeach", "mkanai/WRShd"))
```

## Usage

1. First, you need to prepare the 1000 Genomes Project Phase 3 (version 5) datasets. Our provided script `init_data.sh` will retreive the vcf files from their ftp site and convert them to appropriate PLINK bed files. If you've already downloaded the vcf files, add `--skip-download` option to skip downloading.

2. Second, run a simulation. The script `run.sh` will conduct `$start - $end` permutations and measure the minimum *P*-values of each ancestry/meta-analysis for every permutation. The result will be saved under `analysis/$MAF/ALL` directory. We conducted 100,000 permutations for our main analysis (Kanai, M. *et al.*). but please note that **this step requires huge computational time.**

3. After running all permutations, you can compute and plot the empirical significance thresholds with `utils/plot.R`.

```{bash}
# fetch the 1000 Genomes Phase3 data and convert them appropriately.
./init_data.sh

# run a simulation.
# ./run.sh $maf $start $end $N
# e.g. run 100,000 permutaions for the MAF005 dataset, by splitting the jobs into 100 sub-units.
./run.sh 005 1 100000 100

# plot
# Rscript utils/plot.R $maf $start $end
Rscript utils/plot.R MAF005 1 100000
```

## Citation
* Kanai, M. *et al.* Empirical estimation of genome-wide significance thresholds based on the 1000 Genomes Project. (under revision).

