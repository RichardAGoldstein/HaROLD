package rag.cluster;

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.File;

@Command(name = "richards-haplotype-model", footer = "Copyright (c) 2018 Richard A Goldstein", description="", version = "1.0")
public class Options {
    @Option(names = {"-v", "--verbose"}, description = "")
    boolean verbose = false;

    @Option(names = {"-c", "--count-file"}, required = true, description = "")
    File countFile;

    @Option(names = {"-n", "--haplotypes"}, required = true, description = "")
    int haplotypes;

    @Option(names = {"-g", "--gamma-cache"}, description = "", hidden = true)
    int gammaCache = 0;

    @Option(names = { "-V", "--version" }, versionHelp = true, description = "")
    boolean versionRequested;

    @Option(names = { "-h", "-?", "--help" }, usageHelp = true, description = "give this help list")
    protected boolean helpRequested;

/*    @CommandLine.Parameters(arity = "1..*", paramLabel = "FILE", description = "one ore more files")
    File[] files = null;
    */
}

