package rag.cluster;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import picocli.CommandLine;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

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
                System.out.printf("Main: arguments = %s\n", String.join(" ", args));

                validateOptions(options);

                // Setup
                System.out.printf("Main: seed = %d\n", options.randomSeed);
                GammaCalc gammaCalc = GammaCalc.get(options.gammaCache);

                // fraction of sites to use when optimising alpha parameters
                Constants.USE_FRAC[0] = options.alpha_frac;
                Constants.USE_FRAC[1] = options.alpha_frac;

                long fileSeed = options.randomSeed;

                List<Cluster> clusters = new ArrayList<>();
                for (int i = 0; i < options.countFile.length; i++) {
                    Cluster cluster = new Cluster(options.countFile[i],
                            options.haplotypes[i],
                            options.initialAlphaParams,
                            gammaCalc,
                            fileSeed++,
                            options.verbose);
                    cluster.initialise();
                    clusters.add(cluster);
                }

                // Optimise
                optimise(clusters, options);

                long endTime = System.currentTimeMillis();
                System.out.printf("Main: Execution time = %.2fs\n", (endTime - startTime) / 1000.0);
            }

        } catch (CommandLine.ParameterException ex) {
            System.err.println(ex.getMessage());
            ex.getCommandLine().usage(System.err);
        } catch (Exception ex) {
            throw new CommandLine.ExecutionException(cmd, "Error", ex);
        }

    }

    private void validateOptions(Options options) {
        if (options.countFile.length != options.haplotypes.length) {
            String msg = String.format("You have %d files but %d haplotype numbers.\n", options.countFile.length, options.haplotypes.length);
            throw new RuntimeException(msg);
        }
    }

    private void optimise(List<Cluster> clusters, Options options) {
        final ExecutorService threadPool = Executors.newFixedThreadPool(options.threads);

        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, options.tol);

        PointValuePair previous = new PointValuePair(null, Double.NEGATIVE_INFINITY);
        int iteration = 0;

        final double[] currentAlphaParams = options.initialAlphaParams;

        // optimise until convergence
        while (true) {
            iteration++;
            System.out.println("Main: Optimise haplotype frequencies");
            List<Future<Double>> futures = new ArrayList<>();
            for (final Cluster cluster : clusters) {
                // optimise each cluster haplotypes independently (no synchronisation req)
                Future<Double> future = threadPool.submit(cluster::run);
                futures.add(future);
            }
            List<Double> output = Main.getFutureResults(futures);
            double total = output.stream().mapToDouble(Double::doubleValue).sum();
            System.out.printf("Main: Optimised haplotype frequencies; total = %.7f\n", total);

            // optimise the error alpha parameter
            System.out.printf("Main: Optimise alpha; start = [%.3f, %.3f]\n", currentAlphaParams[0], currentAlphaParams[1]);
            double[] tempAlpha = optimiseAlpha(clusters, currentAlphaParams, threadPool);
            currentAlphaParams[0] = tempAlpha[0];
            currentAlphaParams[1] = tempAlpha[1];

            // calculate current lnl
            futures = new ArrayList<>();
            for (final Cluster cluster : clusters) {
                // optimise each cluster haplotypes independently (no synchronisation req)
                Future<Double> future = threadPool.submit(() -> cluster.calculateCurrent(currentAlphaParams));
                futures.add(future);
            }
            output = Main.getFutureResults(futures);
            total = output.stream().mapToDouble(Double::doubleValue).sum();

            System.out.printf("Main: Optimised alpha; [%.3f, %.3f]; total = %.7f\n", currentAlphaParams[0], currentAlphaParams[1], total);

            PointValuePair current = new PointValuePair(null, total);

            if (convergenceChecker.converged(iteration, previous, current)) {
                break;
            }

            previous = current;
        }

        threadPool.shutdown();

        System.out.println("\nMain: Converged.");
        System.out.println("\n\n========================= RESULTS =========================");

        double finalLnl = 0;
        for (Cluster cluster : clusters) {
            System.out.println();
            finalLnl += cluster.printResults();
        }

        System.out.printf("\nMain: Final total likelihood = %.7f\n", finalLnl);
    }

    private double[] optimiseAlpha(List<Cluster> clusters, double[] startAlpha, ExecutorService threadPool) {
        MultivariateFunction clusterAlphaOptimise = new OptimiseAlphaFunction(clusters, threadPool);

        double[] lb_alpha = new double[]{0.1, 0.0001};
        double[] ub_alpha = new double[]{1000.0, 10.0};

        MultivariateOptimizer optimize = new BOBYQAOptimizer(2*2,0.01,1e-6);
        OptimizationData[] optimizationData = new OptimizationData[]{
                new InitialGuess(startAlpha),
                new MaxEval(1000000),
                GoalType.MINIMIZE,
                new ObjectiveFunction(clusterAlphaOptimise),
                new SimpleBounds(lb_alpha, ub_alpha)};

        return optimize.optimize(optimizationData).getPoint();
    }

    private class OptimiseAlphaFunction implements MultivariateFunction {
        final List<Cluster> clusters;
        final ExecutorService threadPool;
        private OptimiseAlphaFunction(final List<Cluster> clusters, ExecutorService threadPool) {
            this.clusters = clusters;
            this.threadPool = threadPool;
        }

        @Override
        public double value(double[] point) {
            List<Future<Double>> futures = new ArrayList<>();
            for (final Cluster cluster : clusters) {
                // optimise each cluster haplotypes independently (no synchronisation req)
                Future<Double> future = threadPool.submit(() -> cluster.optimiseAlpha(1, point));
                futures.add(future);
            }
            List<Double> output = Main.getFutureResults(futures);
            return output.stream().mapToDouble(Double::doubleValue).sum();
        }
    }

    private static <T> List<T> getFutureResults(List<Future<T>> futures) {
        List<T> results = new ArrayList<>();

        for (Future<T> f : futures) {
            try {
                results.add(f.get());
            } catch (Exception e) {
                e.printStackTrace();
                throw new RuntimeException(e);
            }
        }
        return results;
    }

}
