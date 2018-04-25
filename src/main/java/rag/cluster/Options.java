package rag.cluster;

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.File;
import java.util.Arrays;

@Command(name = "richards-haplotype-model", footer = "Copyright (c) 2018 Richard A Goldstein", description = "", version = "1.0")
public class Options {
    @Option(names = {"-v", "--verbose"}, description = "")
    boolean verbose = false;

    // TODO: -c and -n should be n-arity inputs - make arrays
    @Option(names = {"-c", "--count-file"}, arity = "1..*", required = true, description = "file containing list of count files")
    File[] countFile;

    @Option(names = {"-n", "--haplotypes"}, arity = "1..*", required = true, description = "number of haplotypes")
    int[] haplotypes;

    @Option(names = {"-g", "--gamma-cache"}, description = "", hidden = true)
    int gammaCache = 0;

    @Option(names = {"-V", "--version"}, versionHelp = true, description = "")
    boolean versionRequested;

    @Option(names = {"--alpha-frac"}, required = false, description = "Fraction of sites to use to optimise error parameters")
    double alpha_frac = 1.0;

    @Option(names = {"-h", "-?", "--help"}, usageHelp = true, description = "give this help list")
    protected boolean helpRequested;

    @Option(names = {"-s", "--seed"}, description = "")
    long randomSeed = System.currentTimeMillis();

    // TODO: how to set default
    @Option(names = {"-a", "--initial-alpha"}, arity = "2", description = "")
    double[] initialAlphaParams;

    @Option(names = {"--fix-alpha"}, description="Do not optimise the alpha error parameters")
    boolean fixAlpha;

    @Option(names = {"--threads"})
    int threads = 1;

    @Option(names = {"--tol"})
    double tol = Constants.DEFAULT_TOL;
}

