# HaROLD: Haplotype RecOstruction using Longitudinal Data

This program performs haplotype reconstruction on longitudinal deep sequencing samples, by analysing co-varying variants in a probabilistic framework. 

## Installation

HaROLD requires Java 8 (or newer). Download pre-built binaries from [releases](https://github.com/RichardAGoldstein/HaROLD/releases).

Alternatively, you can clone the source and build using the following commands (requires [Maven](https://maven.apache.org/)):

```
git clone https://github.com/RichardAGoldstein/HaROLD.git
cd HaROLD
mvn compile assembly:single
```

This will create a Java JAR file in the 'target' directory.

## Usage

View program options:

```
$ java -jar harold-1.0.jar -h

Usage:

HaROLD haplotype reconstruction program

java -jar harold-1.0.jar [-hvV] [--alpha-frac=<alpha_frac>]
                         [--error-opt-iter=<errorOptimiseIterations>]
                         [--threads=<threads>] [--tol=<tol>] [-g=<gammaCache>]
                         [-s=<randomSeed>] [-a=<initialAlphaParams>
                         <initialAlphaParams>]... -c=<countFile>...
                         [-c=<countFile>...]... -n=<haplotypes>...
                         [-n=<haplotypes>...]...

Description:

HaROLD (HAplotype Reconstruction Of Longitudinal Deep sequencing Data) performs
haplotype reconstruction on longitudinal deep sequencing samples, by analysing
co-varying variants in a probabilistic framework.

HaROLD reads in a set of files, one for each timepoint. These files should be in the
output format of bam-readcounts.

Run using: java -jar harold-1.0.jar -c <count file> -n <no. of haplotypes>

Options:
  -c, --count-file=<countFile>...
                            File containing list of count files
  -n, --haplotypes=<haplotypes>...
                            Number of haplotypes
  -g, --gamma-cache=<gammaCache>
                            Number of Gamma function calculations to cache
  -s, --seed=<randomSeed>   Seed for random number generator
      --threads=<threads>   Number of processors for multi-threaded operation
      --alpha-frac=<alpha_frac>
                            Fraction of sites to use to optimise error parameters
  -a, --initial-alpha=<initialAlphaParams> <initialAlphaParams>
                            Initial parameter values for error model
      --error-opt-iter=<errorOptimiseIterations>
                            Limit error parameter optimisation to n rounds (0 means
                              no limit)
      --tol=<tol>           Optimisation tolerance
  -h, -?, --help            Show this help
  -v, --verbose
  -V, --version             Show version

Example contains a simple example: three different timepoints of synthetically created data, 
consisting of three different mixtures of the two CMV sequences KP745665.1 and KP745692.1.
(See Supplementary Material in bioRxiv document for details.)

e.g.(in the example folder)  java -jar harold-1.0.jar -c filelist -n 2

Program description and documentation available at:

R. A. Goldstein, A. U. Tamuri, S. Roy, J. Breuer, 2018, Haplotype assignment of virus NGS data using co-variation of variant frequencies, bioRxiv doi: 10.1101/444877.


Copyright (c) 2018 Richard A Goldstein
```


