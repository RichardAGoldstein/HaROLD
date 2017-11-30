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
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;


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
                alphaParams[iTimePoint][iHaplo] = (1.0 + 0.1 * (random.nextDouble() - 0.5)) /(nHaplo-iHaplo);
            }
        }
        return alphaParams;       
    }
    
    
    void run() {
        double initialAlpha = 0.001;
        double initialBeta = 1.0;
        double initialS = 0.0001;
        double[][] currentAlphaParams = initialiseAlphaParams();
        double[] currentErrorParams = new double[3];
        currentErrorParams[0] = initialAlpha;
        currentErrorParams[1] = initialBeta;
        currentErrorParams[2] = initialS;
        double trustRadius = 0.01;
        double[] optPoint = null;
        
        for (int iIter = 0; iIter < 10; iIter++) {
            
            dataSet.assignHaplotypes(currentAlphaParams, currentErrorParams);
            
            if (nHaplo > 2) {
                trustRadius = Math.max(1.0E-6, Math.pow(0.1, iIter+2));
                MultivariateOptimizer optimize = new BOBYQAOptimizer(2*nHaplo-2,0.01,trustRadius);
                double[] lb_alpha = new double[nHaplo-1];
                Arrays.fill(lb_alpha, 1.0E-8);
                double[] ub_alpha = new double[nHaplo-1];
                Arrays.fill(ub_alpha, 1.0);
                optPoint = new double[nHaplo-1];
                trustRadius = 1.0E-6;

                for (int iTimePoint = 0; iTimePoint < dataSet.nTimePoints; iTimePoint++) { 
                    System.out.println(Arrays.toString(currentAlphaParams[iTimePoint]) + "\t" + Arrays.toString(lb_alpha) + "\t" + Arrays.toString(ub_alpha));
                    dataSet.setOptType(1, iTimePoint);
                    OptimizationData[] parm = new OptimizationData[]{
                        new InitialGuess(currentAlphaParams[iTimePoint]),
                        new MaxEval(1000000),
                        GoalType.MINIMIZE,
                        new ObjectiveFunction(dataSet),
                        new SimpleBounds(lb_alpha,ub_alpha)};
                    optPoint = optimize.optimize(parm).getPoint();  // It will be 'THE BEST'
                    System.out.println(Arrays.toString(optPoint));
                    for (int iHaplo = 0; iHaplo < nHaplo-1; iHaplo++) {
                        currentAlphaParams[iTimePoint][iHaplo] = optPoint[iHaplo]; 
                    }

                }
            } else {
                double optSinglePoint = 0.0;

                for (int iTimePoint = 0; iTimePoint < dataSet.nTimePoints; iTimePoint++) { 
                    dataSet.setOptType(1, iTimePoint);
                    optSinglePoint = currentAlphaParams[iTimePoint][0];
                    optSinglePoint = fmin(1.0E-8, 1.0, 1.0E-6);
                    System.out.println(optSinglePoint);
                    currentAlphaParams[iTimePoint][0] = optSinglePoint; 
                }
            }
            
            trustRadius = Math.max(1.0E-6, Math.pow(0.1, iIter+2));
            double[] lb_alpha = {1.0E-8, 1.0E-4, 1.0E-8};
            double[] ub_alpha = {0.1, 5.0, 1.0E-4};
            trustRadius = 1.0E-6;
            MultivariateOptimizer optimize = new BOBYQAOptimizer(2*3,0.01,trustRadius);
            dataSet.setOptType(0, 0);
            optPoint = new double[3];
            
            OptimizationData[] parm = new OptimizationData[]{
                new InitialGuess(currentErrorParams),
                new MaxEval(1000000),
                GoalType.MINIMIZE,
                new ObjectiveFunction(dataSet),
                new SimpleBounds(lb_alpha,ub_alpha)};
            optPoint = optimize.optimize(parm).getPoint();  // It will be 'THE BEST'
            System.out.println(Arrays.toString(optPoint));
            currentErrorParams[0] = optPoint[0];
            currentErrorParams[1] = optPoint[1];
            currentErrorParams[2] = optPoint[2];
        }
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
