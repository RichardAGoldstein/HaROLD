package rag.cluster;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;
import picocli.CommandLine;
import sun.jvm.hotspot.types.PointerType;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
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

                // Setup
                System.out.printf("Main: seed - %d\n", options.randomSeed);
                GammaCalc gammaCalc = GammaCalc.get(options.gammaCache);

                List<Cluster> clusters = new ArrayList<>();
                for (int i = 0; i < options.countFile.length; i++) {
                    Cluster cluster = new Cluster(options.countFile[i].getName(),
                            options.haplotypes[i],
                            options.initialAlphaParams,
                            gammaCalc,
                            options.randomSeed,
                            options.verbose);
                    cluster.initialise();
                    clusters.add(cluster);
                }

                // Optimise
                optimise(clusters, options);

                long endTime = System.currentTimeMillis();
                System.out.printf("\nMain: Execution time = %fs", (endTime - startTime) / 1000.0);
            }

        } catch (CommandLine.ParameterException ex) {
            System.err.println(ex.getMessage());
            ex.getCommandLine().usage(System.err);
        } catch (Exception ex) {
            throw new CommandLine.ExecutionException(cmd, "Error", ex);
        }

    }

    private void optimise(List<Cluster> clusters, Options options) {
        final ExecutorService threadPool = Executors.newFixedThreadPool(options.threads);

        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, 1e-3);

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
            double[] tempAlpha = optimiseAlpha(clusters, options, currentAlphaParams, threadPool);
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

        System.out.println("Main: Converged.");
        System.out.println("========================= RESULTS =========================");

        for (Cluster cluster : clusters) {
            System.out.println();
            cluster.printResults();
        }
    }

    private double[] optimiseAlpha(List<Cluster> clusters, Options options, double[] currentAlphaParams, ExecutorService threadPool) {
        MultivariateFunction clusterAlphaOptimise = new OptimiseAlphaFunction(clusters, threadPool);

        double[] lb_alpha = new double[]{0.1, 0.0001};
        double[] ub_alpha = new double[]{1000.0, 10.0};

        MultivariateOptimizer optimize = new BOBYQAOptimizer(2*2,0.01,1e-6);
        OptimizationData[] parm = new OptimizationData[]{
                new InitialGuess(currentAlphaParams),
                new MaxEval(1000000),
                GoalType.MINIMIZE,
                new ObjectiveFunction(clusterAlphaOptimise),
                new SimpleBounds(lb_alpha, ub_alpha)};
        double[] optPoint = optimize.optimize(parm).getPoint();  // Optimise alpha0, alphaE, f

        currentAlphaParams[0] = optPoint[0];  // Update parameters
        currentAlphaParams[1] = optPoint[1];
        return currentAlphaParams;
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
                Future<Double> future = threadPool.submit(new Callable<Double>() {
                    @Override
                    public Double call() throws Exception {
                        return cluster.optimiseAlpha(1, point);
                    }
                });
                futures.add(future);
            }
            List<Double> output = Main.getFutureResults(futures);
            double total = output.stream().mapToDouble(Double::doubleValue).sum();
            return total;
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
