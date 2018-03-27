package rag.cluster;

import picocli.CommandLine;

public class Main {
    public static void main(String[] args) {
        Main m = new Main();
        m.run(args);
    }

    private void run(String[] args) {
        Options options = new Options();
        CommandLine cmd = new CommandLine(options);

        try {
            cmd.parse(args);
            if (cmd.isUsageHelpRequested()) {
                cmd.usage(System.err);
            } else if (cmd.isVersionHelpRequested()) {
                cmd.printVersionHelp(System.err);
            } else {
                long startTime = System.currentTimeMillis();

                System.out.printf("seed: %d\n", options.randomSeed);
                GammaCalc gammaCalc = GammaCalc.get(options.gammaCache);

                if (options.verbose) Cluster.verbose = true;

                Cluster cluster = new Cluster(options.countFile.getName(), options.haplotypes, options.initialAlphaParams, gammaCalc, options.randomSeed, options.verbose);
                cluster.run();

                long endTime = System.currentTimeMillis();

                if (options.verbose) {
                    System.out.printf("Execution time: %fs", (endTime - startTime) / 1000.0);
                }
            }

        } catch (CommandLine.ParameterException ex) {
            System.err.println(ex.getMessage());
            ex.getCommandLine().usage(System.err);
        } catch (Exception ex) {
            throw new CommandLine.ExecutionException(cmd, "Error", ex);
        }

    }
}
