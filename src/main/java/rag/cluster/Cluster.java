package rag.cluster;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;

/**
 *
 * @author rgoldst
 */
public class Cluster {

    private final int nHaplo; // Number of haplotypes, revised based on command line argument
    private int nTimePoints = 0;  // Number of timepoints, revised based on data
    
    private ArrayList<Assignment> assignmentVector = new ArrayList<>();  // Vector of all possible assignments
    private int[] nAssignDiffBases = new int[5]; // Number of assignments with a given number of bases
    private DataSet dataSet;  // Class for holding and manipulating sequence data

    private Random random;
    private boolean verbose; // Print lots of intermediate results

    private int maxIter = 10; // Maximum rounds of optimisation
    private double minImprovement = 1.0;  // Minimum improvement necessary to continue optimisation
    private boolean optimiseAlpha = true;
    private double[] initialAlphaParams;
    
    private double finalLogLikelihood = 0.0;

    private final String name;

    /**
    * Reads in data and initialises
    */  
    Cluster(File countFilesFile, int nHaplo, double[] initialAlpha, GammaCalc gammaCalc, long randomSeed, boolean verbose) {

        this.name = countFilesFile.getName();
        System.out.println(this.name + ": " + countFilesFile.getAbsolutePath());
        this.random = new Random(randomSeed);
        this.verbose = verbose;

        this.initialAlphaParams = initialAlpha;

        // If we've changed the defaults, do not optimise alpha
        if (initialAlpha[0] != Constants.DEFAULT_ALPHA_0 || initialAlpha[1] != Constants.DEFAULT_ALPHA_1) {
            optimiseAlpha = false;
        }

        this.optimiseAlpha = false;

        this.nHaplo = nHaplo;  // Update number of haplotypes
        System.out.printf("%s: haplotypes = %d\n", this.name, this.nHaplo);

        constructAssignments(gammaCalc);  // Construct possible assignments of bases to haplotypes
        dataSet = new DataSet(countFilesFile, nHaplo, assignmentVector, nAssignDiffBases, gammaCalc, random, verbose); // Construct dataset
        nTimePoints = dataSet.getNTimePoints();  // Number of time points in dataset
        System.out.printf("%s: timepoints = %d\n", this.name, this.nTimePoints);
        System.out.printf("%s: sites = %d\n", this.name, dataSet.getSiteCount());
    }

    private double[][] currentHapParams;
    private double[] currentAlphaParams;

    void initialise() {
        this.currentHapParams = initialiseHapParams();  // Start with initial nearly equal haplotype frequencies
        this.currentAlphaParams = Arrays.copyOf(initialAlphaParams, 2);   // Initial values for alpha parameters alpha0 and alphaE
    }

    /**
    * Find best assignments and haplotype frequencies
    */   
    double run() {
        int iIter = 0;

        // System.out.println("Optimising haplotype frequencies");
        // Optimise haplotype frequencies first
        dataSet.updateAllParams(currentHapParams, currentAlphaParams);
        double step1_current_lnl = dataSet.assignHaplotypes();
        dataSet.updateFracConserved();
        double step1_previous_lnl = Double.NEGATIVE_INFINITY;

        while (true) {
            if (Math.abs(step1_current_lnl - step1_previous_lnl) < minImprovement) {
                break;
            }

            if (nHaplo == 2) {   // Simple single parameter optimisation for each time point
                double optSinglePoint;      //  Hapltype frequency parameter
                for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                    dataSet.setOptType(1, iTimePoint, currentHapParams, currentAlphaParams, iIter);    // Tell dataSet what timePoint is being optimised
                    optSinglePoint = fmin(1.0E-8, 1.0, 1.0E-6);    // Find best value within range and tolerance
                    if (verbose) {
                        System.out.println("Optimum piParams\t" + iTimePoint + "\t" + optSinglePoint);  // Output optimum
                    }
                    currentHapParams[iTimePoint][0] = optSinglePoint; // Update current HapParameters
                }
            } else if (nHaplo > 2) {   // Multidimensional parameter optimisation for each time point
                MultivariateOptimizer optimize = new BOBYQAOptimizer(2 * nHaplo - 2, 0.01, 1.0E-6);
                double[] lb_alpha= new double[nHaplo - 1];
                Arrays.fill(lb_alpha, 1.0E-8);  // Lower bound
                double[] ub_alpha = new double[nHaplo - 1];
                Arrays.fill(ub_alpha, 1.0);     // Upper bound

                for (int iTimePoint = 0; iTimePoint < dataSet.nTimePoints; iTimePoint++) {
                    dataSet.setOptType(1, iTimePoint, currentHapParams, currentAlphaParams, iIter);   // Tell dataSet what timePoint is being optimised
                    OptimizationData[] parm = new OptimizationData[]{       // Set up optimisation data
                            new InitialGuess(currentHapParams[iTimePoint]),
                            new MaxEval(1000000),
                            GoalType.MINIMIZE,
                            new ObjectiveFunction(dataSet),
                            new SimpleBounds(lb_alpha, ub_alpha)};
                    double[] optPoint = optimize.optimize(parm).getPoint();  // Optimise
                    if (verbose) {
                        System.out.println("Optimum piParams\t" + iTimePoint + "\t" + Arrays.toString(optPoint));
                    }
                    for (int iHaplo = 0; iHaplo < nHaplo - 1; iHaplo++) {

                        currentHapParams[iTimePoint][iHaplo] = optPoint[iHaplo]; // Update current parameters
                    }
                }
            }

            dataSet.updateAllParams(currentHapParams, currentAlphaParams);
            step1_previous_lnl = step1_current_lnl;
            step1_current_lnl = dataSet.assignHaplotypes();  // Find best set of assignments
        }

        System.out.printf("%s: haplotype frequencies lnl = %.5f\n", this.name, step1_current_lnl);
        return step1_current_lnl;
    }

    double calculateCurrent(double[] currentAlphaParams) {
        this.currentAlphaParams = currentAlphaParams;
        dataSet.updateAllParams(currentHapParams, this.currentAlphaParams);
        return this.dataSet.assignHaplotypes();
    }

    double printResults() {
        dataSet.setOptType(2, 0, currentHapParams, currentAlphaParams, 0);
        dataSet.updateAllParams(currentHapParams, currentAlphaParams);
        finalLogLikelihood = dataSet.assignHaplotypes();  // Find best set of assignments and calculate loglikelihood
        System.out.printf("-------------------- %s --------------------\n", this.name);
        dataSet.printResults();
        return finalLogLikelihood;
    }

    double optimiseAlpha(int iIter, double[] alphaParams) {
        this.currentAlphaParams = alphaParams;
        this.dataSet.setOptType(0, 0, currentHapParams, currentAlphaParams, iIter);   // Instruct dataSet to optimise alpha0 and alphaE
        double val = this.dataSet.value(alphaParams);
        // System.out.printf("optimising alpha: [%.3f, %.3f] -> %.5f\n", alphaParams[0], alphaParams[1], val);
        return val;
    }

    /**
    * Constructs vector of all possible assignments
    */       
    private void constructAssignments(GammaCalc gammaCalc) {
        int nAssignments = pow(Constants.MAX_BASES, nHaplo);  // Theoretical exhaustive number of possible assignments
        for (int iAssign = 0; iAssign < nAssignments; iAssign++) {  // Loop over all possible assignments
            Assignment newAssignment = new Assignment(iAssign, nHaplo, gammaCalc, verbose);
            assignmentVector.add(newAssignment);
            nAssignDiffBases[newAssignment.nPresent]++;
        }
        System.out.printf("%s: assignments = %d\n", name, assignmentVector.size());
    }

    /**
    * Find initial values of parameters representing piHap frequencies of haplotypes
    */
    private double[][] initialiseHapParams() {
        double remaining = 1.0;
        double[][] hapParams = new double[nTimePoints][nHaplo-1];  // Parameters encoding haplotype frequencies
        for (int iHaplo = 0; iHaplo < nHaplo-1; iHaplo++) {
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                hapParams[iTimePoint][iHaplo] = (1.0 + 0.1 * (random.nextDouble() - 0.5)) /(nHaplo-iHaplo); 
            }
        }
        return hapParams;       
    }
    
    
    /**
    * Computes a^b
    */
    private int pow (int a, int b) {  // Computes powers
        if ( b == 0)     return 1;
        if ( b == 1)     return a;
        if (b%2 == 0)    return     pow ( a * a, b/2); //even a=(a^2)^b/2
        else             return a * pow ( a * a, (b-1)/2); //odd  a=a*(a^2)^b/2
    }
    
    /**
    * Brent algorithm for maximising a function in one dimension
    * Adapted from Apache Commons
    */    
    private double fmin (double a, double b, double tol) {
        double c,d,e,eps,xm,p,q,r,tol1,t2, u,v,w,fu,fv,fw,fx,x,tol3;

        c = .5*(3.0 - Math.sqrt(5.0));
        d = 0.0;
        eps = 1.2e-16;
        tol1 = eps + 1.0;
        eps = Math.sqrt(eps);

        v = a + c*(b-a);
        w = v;
        x = v;
        e = 0.0;
        fx=dataSet.value(x);
        fv = fx;
        fw = fx;
        tol3 = tol/3.0;

        xm = .5*(a + b);
        tol1 = eps*Math.abs(x) + tol3;
        t2 = 2.0*tol1;

        while (Math.abs(x-xm) > (t2 - .5*(b-a))) {
            p = q = r = 0.0;
            if (Math.abs(e) > tol1) {
                r = (x-w)*(fx-fv);
                q = (x-v)*(fx-fw);
                p = (x-v)*q - (x-w)*r;
                q = 2.0*(q-r);
                if (q > 0.0) {
                        p = -p;
                } else {
                        q = -q;
                }
                r = e;
                e = d;
            }

            if ((Math.abs(p) < Math.abs(.5*q*r)) && (p > q*(a-x)) && (p < q*(b-x))) {
                d = p/q;
                u = x+d;
                if (((u-a) < t2) || ((b-u) < t2)) {
                    d = tol1;
                    if (x >= xm) d = -d;
                }
            } else {
                if (x < xm) {
                    e = b-x;
                } else {
                    e = a-x;
                }
                d = c*e;
            }

            if (Math.abs(d) >= tol1) {
                u = x+d;
            } else {
                if (d > 0.0) {
                    u = x + tol1;
                } else {
                    u = x - tol1;
                }
            }
            fu = dataSet.value(u);

            if (fx <= fu) {
                if (u < x) {
                    a = u;
                } else {
                    b = u;
                }
            }
            if (fu <= fx) {
                if (u < x) {
                    b = x;
                } else {
                    a = x;
                }
                v = w;
                fv = fw;
                w = x;
                fw = fx;
                x = u;
                fx = fu;
                xm = .5*(a + b);
                tol1 = eps*Math.abs(x) + tol3;
                t2 = 2.0*tol1;
            } else {
                if ((fu <= fw) || (w == x)) {
                    v = w;
                    fv = fw;
                    w = u;
                    fw = fu;
                    xm = .5*(a + b);
                    tol1 = eps*Math.abs(x) + tol3;
                    t2 = 2.0*tol1;
                } else if ((fu > fv) && (v != x) && (v != w)) {
                    xm = .5*(a + b);
                    tol1 = eps*Math.abs(x) + tol3;
                    t2 = 2.0*tol1;
                } else {
                    v = u;
                    fv = fu;
                    xm = .5*(a + b);
                    tol1 = eps*Math.abs(x) + tol3;
                    t2 = 2.0*tol1;
                }
            }
        }
        return x;
    }

}
