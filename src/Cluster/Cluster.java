package Cluster;

import java.util.Arrays;
import java.util.Random;
import java.util.Vector;
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
    int nHaplo = 3; // Number of haplotypes, revised based on command line arguement
    int nTimePoints = 0;  // Number of timepoints, revised based on data
    
    Vector<Assignment> assignmentVector = new Vector<>();  // Vector of all possible assignments
    DataSet dataSet = null;  // Class for holding and manipulating sequence data
    static Random random = new Random(435027);  
    
    static boolean verbose = false; // Print lots of intermediate results
    static double[] useFrac = {0.01, 0.1};  // What fraction of sites to use for global parameters (chosen randomly) 
                                            // First number is for first iteration, second is for later iterations
    
    int maxIter = 10; // Maximium rounds of optimisation
    double minImprovement = 1.0;  // Minimum improvement necessary to continue optimisation
    boolean optimiseAlpha = true;
    double[] initialAlphaParams = {100.0, 0.2, 0.9};
    
    double finalLogLikelihood = 0.0;

    /**
     * @param args File containing list of files and number of haplotypes
     */
    public static void main(String[] args) {
        Cluster cluster = new Cluster(args);
        cluster.run();
    }

    /**
    * Reads in data and initialises
    */  
    Cluster(String[] args) {
        if (args.length < 2) {
            System.out.println("First argument is file containing list of count files");
            System.out.println("Second argument is number of haplotypes");
            System.exit(1);
        }
        if (args.length == 5) {
            initialAlphaParams[0] = Double.parseDouble(args[2]);
            initialAlphaParams[1] = Double.parseDouble(args[3]);
            initialAlphaParams[2] = Double.parseDouble(args[4]);
            optimiseAlpha = false;
        }
        nHaplo = Integer.parseInt(args[1]);  // Update number of haplotypes
        constructAssignments();  // Construct possible assignments of bases to haplotypes
        dataSet = new DataSet(args[0], nHaplo, assignmentVector, useFrac); // Construct dataset
        nTimePoints = dataSet.getNTimePoints();  // Number of time points in dataset
    }
    

    
    /**
    * Find best assignments and haplotype frequencies
    */   
    void run() {
        double[][] currentHapParams = initialiseHapParams();  // Start with initial nearly equal haplotype frequencies
        double[] currentAlphaParams = Arrays.copyOf(initialAlphaParams, 3);   // Initial values for alpha parameters alpha0 and alphaE
        double[] optPoint = null;    // Array for parameters
        double[] lb_alpha = null;    // Array for lower bounds
        double[] ub_alpha = null;    // Array for upper bounds
        MultivariateOptimizer optimize = null;    // Optimiser
        OptimizationData[] parm = null;   // Optimisation parameters
        
        int iIter = 0;
        boolean finished = false;
        double previousLogLikelihood = 1.0E20;
        double currentLogLikelihood = 1.0E20;
        
        while (!finished) {      
            dataSet.updateAllParams(currentHapParams, currentAlphaParams);
            previousLogLikelihood = currentLogLikelihood;
            currentLogLikelihood = dataSet.assignHaplotypes();  // Find best set of assignments
            if (Math.abs(previousLogLikelihood - currentLogLikelihood) < minImprovement) {
                finished = true;
            } else {
                // Optimise haplotype frequencies first           
                if (nHaplo == 2) {   // Simple single parameter optimisation for each time point
                    double optSinglePoint = 0.5;      //  Hapltype frequency parameter
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
                    optPoint = new double[nHaplo-1];

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
                    lb_alpha = new double[3];
                    lb_alpha[0] = 0.1;      // lower bound of alpha0
                    lb_alpha[1] = 0.0001;   // lower bound of alphaE
                    lb_alpha[2] = 0.01;     // lower bound of f
                    ub_alpha = new double[3];  
                    ub_alpha[0] = 1000.0;  // upper bound of alpha0
                    ub_alpha[1] = 10.0;    // upper bound of alphaE
                    ub_alpha[2] = 0.9;     // upper bound of f
                    optimize = new BOBYQAOptimizer(2*3,0.01,1.0E-6);

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
                    currentAlphaParams[2] = optPoint[2];
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
    void constructAssignments() {
        int nAssignments = pow(maxBases, nHaplo);  // Theoretical exhaustive number of possible assignments
        for (int iAssign = 0; iAssign < nAssignments; iAssign++) {  // Loop over all possible assignments
            Assignment newAssignment = new Assignment(iAssign, nHaplo);
            assignmentVector.add(newAssignment);      
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
