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
    static int maxBases = 4; // Maximum number of different bases
    int nHaplo = 3; // Number of haplotypes
//    int[][] assign = null; // different possible assignments of bases to haplotypes
    int nTimePoints = 0;
    
    Vector<Assignment> assignmentVector = new Vector<>();
    DataSet dataSet = null;

    static Random random = new Random(435027);
    static boolean verbose = true; // print lots of intermediate results
    static double useFrac = 1.0;  // what fraction of sites to use (chosen randomly)
    

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Cluster cluster = new Cluster(args);
        cluster.run();
    }

    
    Cluster(String[] args) {
        if (args.length < 2) {
            System.out.println("First argument is file containing list of count files");
            System.out.println("Second argument is number of haplotypes");
            System.exit(1);
        }
        nHaplo = Integer.parseInt(args[1]);  // number of haplotypes
        constructAssignments();  // Construct possible assignments of bases to haplotypes
        dataSet = new DataSet(args[0], nHaplo, assignmentVector); // Construct dataset
        nTimePoints = dataSet.nTimePoints;
    }
    
    double[][] initialiseAlphaParams() {
        double remaining = 1.0;
        double[][] alphaParams = new double[nTimePoints][nHaplo-1];
        for (int iHaplo = 0; iHaplo < nHaplo-1; iHaplo++) {
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                alphaParams[iTimePoint][iHaplo] = (1.0 + 0.1 * (random.nextDouble() - 0.5)) * 1.0/(nHaplo-iHaplo);
            }
        }
        return alphaParams;       
    }
    
    
    void run() {
        double initialEpsilon = 0.1;
        double initialS = 0.0001;
        double[][] currentAlphaParams = initialiseAlphaParams();
        double[] currentErrorParams = new double[2];
        currentErrorParams[0] = initialEpsilon;
        currentErrorParams[1] = initialS;
        double trustRadius = 0.01;
        double[] optPoint = null;
        
        for (int iIter = 0; iIter < 10; iIter++) {
            
            dataSet.assignHaplotypes(currentAlphaParams, currentErrorParams);
//            
//            trustRadius = Math.max(1.0E-6, Math.pow(0.1, iIter+2));
//            MultivariateOptimizer optimize = new BOBYQAOptimizer(2*nHaplo-2,0.01,trustRadius);
//            double[] lb_alpha = new double[nHaplo-1];
//            Arrays.fill(lb_alpha, 1.0E-8);
//            double[] ub_alpha = new double[nHaplo-1];
//            Arrays.fill(ub_alpha, 1.0);
//            optPoint = new double[nHaplo-1];
//            
//            for (int iTimePoint = 0; iTimePoint < dataSet.nTimePoints; iTimePoint++) { 
//                dataSet.setOptType(1, iTimePoint);
//                OptimizationData[] parm = new OptimizationData[]{
//                    new InitialGuess(currentAlphaParams),
//                    new MaxEval(1000000),
//                    GoalType.MINIMIZE,
//                    new ObjectiveFunction(dataSet),
//                    new SimpleBounds(lb_alpha,ub_alpha)};
//                optPoint = optimize.optimize(parm).getPoint();  // It will be 'THE BEST'
//                System.out.println(Arrays.toString(optPoint));
//                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
//                    currentAlphaParams[iHaplo] = optPoint[iHaplo]; 
//                }
//                
            }
//            
//            optimize = new BOBYQAOptimizer(4,0.01,trustRadius);
//            lb_alpha = new double[2];
//            Arrays.fill(lb_alpha, 1.0E-8);
//            ub_alpha = new double[2];
//            Arrays.fill(ub_alpha, 1.0);
//            dataSet.setOptType(0, 0);
//            double[] initial = new double[2];
//            initial[0] = dataSet.currentEpsilon;
//            initial[1] = dataSet.currentS;
//            optPoint = new double[2];
//            
//            OptimizationData[] parm = new OptimizationData[]{
//                new InitialGuess(initial),
//                new MaxEval(1000000),
//                GoalType.MINIMIZE,
//                new ObjectiveFunction(dataSet),
//                new SimpleBounds(lb_alpha,ub_alpha)};
//            optPoint = optimize.optimize(parm).getPoint();  // It will be 'THE BEST'
//            System.out.println(Arrays.toString(optPoint));
//            dataSet.currentEpsilon = optPoint[0];
//            dataSet.currentS = optPoint[1];
//        }
    }

       
    void constructAssignments() {
        // Construct all possible ways of assigning positions in haplotypes to specific bases
        int nAssignments = pow(maxBases, nHaplo);  // Theoretical exhaustive number of possible assignments
        for (int iAssign = 0; iAssign < nAssignments; iAssign++) {  // Loop over all possible assignments
            Assignment newAssignment = new Assignment(iAssign, nHaplo);
            assignmentVector.add(newAssignment);      
        } 
    }
//    
//    void printColHeaders(int iStage, int nTimePoints) {
//        if (iStage == 0) {
//            System.out.println("zzz\nzzz\titer\t  logLike\t     alpha_e       beta_e\t    S   \t    F0");
//        } else if (iStage == 1) {
//            System.out.print("zzz\nzzz\titer\t  logLike\t");
//            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
//                System.out.print("  ");
//                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
//                    System.out.print(" alpha(" + iTimePoint + "," + iHaplo + ")");
//                }
//            }
//            System.out.println();
//        }
//    }
    
    int pow (int a, int b) {  // Computes powers
        if ( b == 0)     return 1;
        if ( b == 1)     return a;
        if (b%2 == 0)    return     pow ( a * a, b/2); //even a=(a^2)^b/2
        else             return a * pow ( a * a, (b-1)/2); //odd  a=a*(a^2)^b/2
    }

}
