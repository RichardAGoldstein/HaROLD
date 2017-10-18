package Cluster;

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


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author rgoldst
 */
public class Cluster {
    
    int nHaplo = 3; // Number of haplotypes
    int nAssignments = 0; // number of possible assignments
    int[][] assign = null; // different possible assignments of bases to haplotypes
    boolean addFlat = false;  // Add a 'garbage' model for random outliers *** NOT IMPLEMENTED***

    double[] lb_glob = {1.0E-8, 1.0E-8, 1.0E-8};  // lower bound on global params (alpha_e, S, F0)
    double[] ub_glob = {0.1, 0.01, 1.0};  // upper bound on global params
    double[] lb_alpha = null;  // lower bound on alphas
    double[] ub_alpha = null;  // upper bound on alphas
    double alpha_e = 0.1;
    double beta_e = 10.0;
    double F0 = 0.02;
    double S = 0.0001;
//    double[][] alpha = null;
    double[][] alpha = {{0.3904,1.889},{4.3321,8.4047},{6.9142,11.7717},{9.8596,12.9925},{6.6161,10.8112},
            {9.5218,10.8302},{13.3333,9.0419},{11.9004,6.5896},{20.9773,4.6513},{22.4221,3.0989}};
    
    static Random random = new Random();
    static boolean verbose = true; // print lots of intermediate results
    static double useFrac = 1.0;  // what fraction of sites to use (chosen randomly)
    int maxRound = 50;  // maximum number of global iterations
    double minDeltaLL = -1.0; // minimum change in log likelihood before stopping
    

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Cluster clus = new Cluster(args);
        clus.run(args);
    }

    Cluster(String[] args) {
        if (args.length == 0) {
            System.out.println("First argument is file containing list of count files");
            System.out.println("First argument is number of haplotypes");
            System.exit(1);
        }
        nHaplo = Integer.parseInt(args[1]);  // number of haplotypes
        constructAssignments();   
    }
    
    
    void run(String[] args) {
        // Initialise stuff
        DataSet dataSet = new DataSet(args[0], nHaplo, nAssignments, assign, addFlat); // Construct dataset
//        alpha = new double[dataSet.getNTimePoints()][nHaplo];
        lb_alpha = new double[nHaplo];
        Arrays.fill(lb_alpha, 1.0E-8);
        ub_alpha = new double[nHaplo];
        Arrays.fill(ub_alpha, 100.0);
//        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
//            for (int iTimePoint = 0; iTimePoint < dataSet.getNTimePoints(); iTimePoint++) {
//                alpha[iTimePoint][iHaplo] = 2.0 * random.nextDouble();
//            }
//        }
        
        
        
        // Optimise parameters
        optimiseModel(dataSet);  
    }
    
    void optimiseModel(DataSet dataSet) {
        printColHeaders(0, 0);
        dataSet.setStage(0, alpha_e, beta_e, S, F0, alpha);  // Adjust global parameters and E step
        double previous = adjustGlobals(dataSet);  // Keep track of improvements
        double current = previous + 1000.0;
        
        int iRound_0 = 0;  // loop over global optimisations
        boolean finished = false;
        while (!finished) {
            printColHeaders(1, dataSet.getNTimePoints());
            //  loop over EM steps for inner loop
            double valueAlpha = 0.0;
            dataSet.setStage(1,alpha_e, beta_e, S, F0, alpha);
            for (int iRound_1 = 0; iRound_1 < 50; iRound_1++) {
                valueAlpha = 0.0;
                // M step for each timepoint
                for (int iTimePoint = 0; iTimePoint < dataSet.getNTimePoints(); iTimePoint++) {
                    dataSet.setITimePoint(iTimePoint);
                    valueAlpha += adjustAlphas(dataSet, iTimePoint);
                }
                // E step
                dataSet.reAdjust(alpha_e, beta_e, S, F0);
                
                if (iRound_1%10 == 0) {
                    System.out.format("zzz\t%d\t%10.4f   ", iRound_1, valueAlpha);
                    for (int iTimePoint = 0; iTimePoint < dataSet.getNTimePoints(); iTimePoint++) {
                        System.out.print("    ");
                        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                            System.out.format("  %8.4f", alpha[iTimePoint][iHaplo]);
                        }
                    }
                    System.out.println();      
                }
            }
            
            System.out.format("zzz\tConv'g\t%10.4f   ", valueAlpha);
            for (int iTimePoint = 0; iTimePoint < dataSet.getNTimePoints(); iTimePoint++) {
                System.out.print("    ");
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    System.out.format("  %8.4f", alpha[iTimePoint][iHaplo]);
                }
            }
            System.out.println();
            
            printColHeaders(0, 0);
            dataSet.setStage(0, alpha_e, beta_e, S, F0, alpha); // Adjust global parameters and E step
            previous = current;
            current = adjustGlobals(dataSet);
            iRound_0++;
            if (iRound_0 > maxRound || (current - previous)>minDeltaLL) {
                finished = true;
            }
        }
        dataSet.printHaplotypes(alpha_e, beta_e, S, F0);
    }

    double adjustAlphas(DataSet dataSet, int iTimePoint) {  // Find best parameters for understanding variability
        MultivariateOptimizer optimize = new BOBYQAOptimizer(2*nHaplo,0.01,1.0E-6);
        OptimizationData[] parm = new OptimizationData[]{
            new InitialGuess(alpha[iTimePoint]),
            new MaxEval(1000000),
            GoalType.MINIMIZE,
            new ObjectiveFunction(dataSet),
            new SimpleBounds(lb_alpha,ub_alpha)};
        alpha[iTimePoint] = optimize.optimize(parm).getPoint();  // It will be 'THE BEST'
        return dataSet.value(alpha[iTimePoint]);
    }
    
    
    double adjustGlobals(DataSet dataSet) {  // Find best parameters for understanding variability
        int nVar = 3;  // set up and initialise parameters
        double[] params = new double[3];
        params[0] = alpha_e;
        params[1] = S;
        params[2] = F0;

        // Optimise 
        MultivariateOptimizer optimize = new BOBYQAOptimizer(2*nVar,0.01,1.0E-6);
        OptimizationData[] parm = new OptimizationData[]{
            new InitialGuess(params),
            new MaxEval(1000000),
            GoalType.MINIMIZE,
            new ObjectiveFunction(dataSet),
            new SimpleBounds(lb_glob,ub_glob)};
        double[] best = optimize.optimize(parm).getPoint();  // It will be 'THE BEST'
        alpha_e = best[0];
        S = best[1];
        F0 = best[2];
        double value = dataSet.value(best);
        System.out.format("zzz\tConv'd\t%10.4f\t%12.6f  %12.6f\t%10.8f\t%10.8f\n",
            value, alpha_e, beta_e, S, F0);
        return value;
    }
       
    void constructAssignments() {
        int[] powersOfTwo = new int[nHaplo+1];
        powersOfTwo[0]=1;
        for (int iPow = 1; iPow < nHaplo+1; iPow++) {
            powersOfTwo[iPow] = powersOfTwo[iPow-1]*2;
        }
        nAssignments = powersOfTwo[nHaplo];  // All ways of distributing variants amongst haplotypes
        assign = new int[nAssignments][nHaplo];    //  List of possible assignments
        for (int iAssignment = 0; iAssignment<nAssignments; iAssignment++) {
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                if (((iAssignment) & powersOfTwo[iHaplo]) > 0) assign[iAssignment][iHaplo] = 1;  // generating all possibilities
            }
        } 
        if (verbose) {
            System.out.println("zzz\tVarious assignments");
            for (int iAssignment = 0; iAssignment < nAssignments; iAssignment++) {
                System.out.println("zzz\t" + iAssignment + "\t" + Arrays.toString(assign[iAssignment]));
            }
            System.out.println("zzz");
        }       
    }
    
    void printColHeaders(int iStage, int nTimePoints) {
        if (iStage == 0) {
            System.out.println("zzz\nzzz\titer\t  logLike\t     alpha_e       beta_e\t    S   \t    F0");
        } else if (iStage == 1) {
            System.out.print("zzz\nzzz\titer\t  logLike\t");
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                System.out.print("  ");
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    System.out.print(" alpha(" + iTimePoint + "," + iHaplo + ")");
                }
            }
            System.out.println();
        }
    }
    
    
}
