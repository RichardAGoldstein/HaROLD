# HaROLD: Haplotype RecOstruction using Longitudinal Data

This program performs haplotype reconstruction on longitudinal deep sequencing samples, by analysing co-varying variants in a probabilistic framework. 

## Installation

You can download the pre-built releases from https://github.com/RichardAGoldstein/HaROLD/releases

HaROLD requires Java 8 (or newer).

Alternatively, you can clone the source and build using Maven:

```
git clone https://github.com/RichardAGoldstein/HaROLD.git
cd HaROLD
mvn
```

This will create a JAR file in the target directory.

## Usage

View program options:

```
$ java -jar target/harold-1.0-SNAPSHOT-jar-with-dependencies.jar -h

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
Copyright (c) 2018 Richard A Goldstein

~/Documents/2018/haplotypes/HaROLD
```

