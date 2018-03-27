package rag.cluster;

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

    static int maxBases = 4; // Maximum number of different bases
    final int nHaplo; // Number of haplotypes, revised based on command line argument
    int nTimePoints = 0;  // Number of timepoints, revised based on data
    
    ArrayList<Assignment> assignmentVector = new ArrayList<>();  // Vector of all possible assignments
    int[] nAssignDiffBases = new int[5]; // Number of assignments with a given number of bases
    DataSet dataSet = null;  // Class for holding and manipulating sequence data
    // static Random random = new Random(435027);

    static Random random;

    static boolean verbose; // Print lots of intermediate results
    static double[] useFrac = {0.01, 0.1};  // What fraction of sites to use for global parameters (chosen randomly) 
                                            // First number is for first iteration, second is for later iterations
    
    int maxIter = 10; // Maximum rounds of optimisation
    double minImprovement = 1.0;  // Minimum improvement necessary to continue optimisation
    boolean optimiseAlpha = true;
    double[] initialAlphaParams;
    
    double finalLogLikelihood = 0.0;


    /**
    * Reads in data and initialises
    */  
    Cluster(String countFilesFile, int nHaplo, double[] initialAlpha, GammaCalc gammaCalc, long randomSeed, boolean verbose) {
        Cluster.random = new Random(randomSeed);
        Cluster.verbose = verbose;

        this.initialAlphaParams = initialAlpha;

        // If we've changed the defaults, do not optimise alpha
        if (initialAlpha[0] != Constants.DEFAULT_ALPHA_0 || initialAlpha[1] != Constants.DEFAULT_ALPHA_1) {
            optimiseAlpha = false;
        }

        this.nHaplo = nHaplo;  // Update number of haplotypes
        System.out.printf("haplotypes: %d\n", this.nHaplo);

        constructAssignments(gammaCalc);  // Construct possible assignments of bases to haplotypes
        dataSet = new DataSet(countFilesFile, nHaplo, assignmentVector, nAssignDiffBases, useFrac, gammaCalc); // Construct dataset
        nTimePoints = dataSet.getNTimePoints();  // Number of time points in dataset
    }
    

    /**
    * Find best assignments and haplotype frequencies
    */   
    void run() {
        double[][] currentHapParams = initialiseHapParams();  // Start with initial nearly equal haplotype frequencies
        double[] currentAlphaParams = Arrays.copyOf(initialAlphaParams, 2);   // Initial values for alpha parameters alpha0 and alphaE
        double[] optPoint;    // Array for parameters
        double[] lb_alpha;    // Array for lower bounds
        double[] ub_alpha;    // Array for upper bounds
        MultivariateOptimizer optimize = null;    // Optimiser
        OptimizationData[] parm;   // Optimisation parameters
        
        int iIter = 0;
        boolean finished = false;
        double previousLogLikelihood = 1.0E20;
        double currentLogLikelihood = 1.0E20;
        
        while (!finished) {      
            dataSet.updateAllParams(currentHapParams, currentAlphaParams);
            previousLogLikelihood = currentLogLikelihood;
            currentLogLikelihood = dataSet.assignHaplotypes();  // Find best set of assignments
            dataSet.updateFracConserved();
            if (Math.abs(previousLogLikelihood - currentLogLikelihood) < minImprovement) {
                finished = true;
            } else {
                // Optimise haplotype frequencies first           
                if (nHaplo == 1) {
                } else if (nHaplo == 2) {   // Simple single parameter optimisation for each time point
                    double optSinglePoint;      //  Hapltype frequency parameter
                    for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) { 
                        dataSet.setOptType(1, iTimePoint, currentHapParams, currentAlphaParams, iIter);    // Tell dataSet what timePoint is being optimised
                        optSinglePoint = fmin(1.0E-8, 1.0, 1.0E-6);    // Find best value within range and tolerance
                        if (verbose) {
                            System.out.println("Optimum piParams\t" + iTimePoint + "\t" + optSinglePoint);  // Output optimum
                        }
                        currentHapParams[iTimePoint][0] = optSinglePoint; // Update current HapParameters
                    }
                } else {   // Multidimensional parameter optimisation for each time point
                    optimize = new BOBYQAOptimizer(2*nHaplo-2,0.01,1.0E-6);
                    lb_alpha = new double[nHaplo-1];
                    Arrays.fill(lb_alpha, 1.0E-8);  // Lower bound
                    ub_alpha = new double[nHaplo-1];
                    Arrays.fill(ub_alpha, 1.0);     // Upper bound

                    for (int iTimePoint = 0; iTimePoint < dataSet.nTimePoints; iTimePoint++) {
                        dataSet.setOptType(1, iTimePoint, currentHapParams, currentAlphaParams, iIter);   // Tell dataSet what timePoint is being optimised
                        parm = new OptimizationData[]{       // Set up optimisation data
                            new InitialGuess(currentHapParams[iTimePoint]),
                            new MaxEval(1000000),
                            GoalType.MINIMIZE,
                            new ObjectiveFunction(dataSet),
                            new SimpleBounds(lb_alpha,ub_alpha)};
                        optPoint = optimize.optimize(parm).getPoint();  // Optimise
                        if (verbose) {
                            System.out.println("Optimum piParams\t" + iTimePoint + "\t" + Arrays.toString(optPoint));
                        }
                        for (int iHaplo = 0; iHaplo < nHaplo-1; iHaplo++) {
                            currentHapParams[iTimePoint][iHaplo] = optPoint[iHaplo]; // Update current parameters
                        }

                    }
                }

                if (optimiseAlpha) {
                    dataSet.setOptType(0, 0, currentHapParams, currentAlphaParams, iIter);   // Instruct dataSet to optimise alpha0 and alphaE   
                    lb_alpha = new double[2];
                    lb_alpha[0] = 0.1;      // lower bound of alpha0
                    lb_alpha[1] = 0.0001;   // lower bound of alphaE
                    ub_alpha = new double[2];  
                    ub_alpha[0] = 1000.0;  // upper bound of alpha0
                    ub_alpha[1] = 10.0;    // upper bound of alphaE
                    optimize = new BOBYQAOptimizer(2*2,0.01,1.0E-6);

                    parm = new OptimizationData[]{
                        new InitialGuess(currentAlphaParams),
                        new MaxEval(1000000),
                        GoalType.MINIMIZE,
                        new ObjectiveFunction(dataSet),
                        new SimpleBounds(lb_alpha,ub_alpha)};
                    optPoint = optimize.optimize(parm).getPoint();  // Optimise alpha0, alphaE, f
                    if (verbose) {
                        System.out.println("Optimimum alphaParams:\t" + Arrays.toString(optPoint));
                    }
                    currentAlphaParams[0] = optPoint[0];  // Update parameters
                    currentAlphaParams[1] = optPoint[1];
                }
                
                iIter++;
                if (iIter == maxIter) {
                    finished = true;
                }
            }
        }        
        dataSet.setOptType(2, 0, currentHapParams, currentAlphaParams, 0);
        dataSet.updateAllParams(currentHapParams, currentAlphaParams);
        finalLogLikelihood = dataSet.assignHaplotypes();  // Find best set of assignments and calculate loglikelihood
        dataSet.printResults();
    }

    /**
    * Constructs vector of all possible assignments
    */       
    void constructAssignments(GammaCalc gammaCalc) {
        int nAssignments = pow(maxBases, nHaplo);  // Theoretical exhaustive number of possible assignments
        for (int iAssign = 0; iAssign < nAssignments; iAssign++) {  // Loop over all possible assignments
            Assignment newAssignment = new Assignment(iAssign, nHaplo, gammaCalc);
            assignmentVector.add(newAssignment);
            nAssignDiffBases[newAssignment.nPresent]++;
        } 
    }
    
    /**
    * Find initial values of parameters representing piHap frequencies of haplotypes
    */
    double[][] initialiseHapParams() { 
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
    int pow (int a, int b) {  // Computes powers
        if ( b == 0)     return 1;
        if ( b == 1)     return a;
        if (b%2 == 0)    return     pow ( a * a, b/2); //even a=(a^2)^b/2
        else             return a * pow ( a * a, (b-1)/2); //odd  a=a*(a^2)^b/2
    }
    
    /**
    * Brent algorithm for maximising a function in one dimension
    * Adapted from Apache Commons
    */    
    double fmin (double a, double b, double tol) {
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
