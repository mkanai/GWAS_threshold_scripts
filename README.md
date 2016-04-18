# Empirical estimation of genome-wide significance thresholds based on the 1000 Genomes Project.

This repository contains the scripts used in Kanai, M. *et al.*

## Requirements
* Python
* R
* [GNU Parallel](http://www.gnu.org/software/parallel/)
* Sun Grid Engine (or equivalent) environment
    * We use `qsub` command in our scripts. For other environments, please replace the command appropriately.
* Our customized version of [PLINK 1.9](https://www.cog-genomics.org/plink2/)
    * We added several features to conduct our GWAS simulations against randomly assigned case-control phenotypes.
    * Our customized version is based on `PLINK v1.90b3l` (chrchang/plink-ng@1d8cbd4).
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

1. First, you need to prepare the 1000 Genomes Project Phase 3 (version 5) datasets. Our provided script `init_data.sh` will retreive and convert their vcf files to appropriate PLINK bed files. If you've already fetched the vcf files, add `--skip-download` option to skip downloading.

2. Second, run a simulation. The script `run.sh` will conduct `$start - $end` permutations and measure the minimum P-values of each ancestry/meta-analysis for every permutation. The result will be saved under `analysis/$MAF/ALL` directory.

3. Finally, ...

```{bash}
./init_data.sh

./run.sh $start $end
```

## Citation
* Kanai, M. *et al.* Empirical estimation of genome-wide significance thresholds based on the 1000 Genomes Project. (under revision).
