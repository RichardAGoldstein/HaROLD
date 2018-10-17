package rag.harold;

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.File;

@Command(name = "HaROLD",
        footer = "Copyright (c) 2018 Richard A Goldstein",
        description = "HaROLD (HAplotype Reconstruction Of Longitudinal Deep sequencing Data) " +
                "performs haplotype reconstruction on longitudinal deep sequencing samples, by " +
                "analysing co-varying variants in a probabilistic framework.",
        version = "1.0",
        sortOptions = false,
        headerHeading = "Usage:%n%n",
        synopsisHeading = "%n",
        descriptionHeading = "%nDescription:%n%n",
        parameterListHeading = "%nParameters:%n",
        optionListHeading = "%nOptions:%n",
        header = "HaROLD haplotype reconstruction")
public class Options {
    @Option(names = {"-c", "--count-file"}, arity = "1..*", required = true, description = "File containing list of count files")
    File[] countFile;

    @Option(names = {"-n", "--haplotypes"}, arity = "1..*", required = true, description = "Number of haplotypes")
    int[] haplotypes;

    @Option(names = {"-g", "--gamma-cache"}, description = "Number of Gamma function calculations to cache")
    int gammaCache = 0;

    @Option(names = {"-s", "--seed"}, description = "Seed for random number generator")
    long randomSeed = System.currentTimeMillis();

    @Option(names = {"--threads"}, description = "Number of processors for multi-threaded operation")
    int threads = 1;

    @Option(names = {"--alpha-frac"}, required = false, description = "Fraction of sites to use to optimise error parameters")
    double alpha_frac = 1.0;

    @Option(names = {"-a", "--initial-alpha"}, arity = "2", description = "Initial parameter values for error model")
    double[] initialAlphaParams = new double[]{Constants.DEFAULT_ALPHA_0, Constants.DEFAULT_ALPHA_1};

    @Option(names = {"--error-opt-iter"}, arity = "1", description = "Limit error parameter optimisation to n rounds (0 means no limit)")
    int errorOptimiseIterations = 0;

    @Option(names = {"--tol"}, description = "Optimisation tolerance")
    double tol = Constants.DEFAULT_TOL;

    @Option(names = {"-h", "-?", "--help"}, usageHelp = true, description = "Show this help")
    protected boolean helpRequested;

    @Option(names = {"-v", "--verbose"}, description = "")
    boolean verbose = false;

    @Option(names = {"-V", "--version"}, versionHelp = true, description = "Show version")
    boolean versionRequested;
}

