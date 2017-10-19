/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Cluster;

import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Random;
import java.util.Vector;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.CombinatoricsUtils;

/**
 * Container for the data at a single site
 * 
 * @author rgoldst
 */
public class Site {
    
    int iSite;
    int nTimePoints;
    int nHaplo = 0;
    int maxBases = 2;
    int nBases = 0;
    int[] nAssignments = null;
    int[][] assign = null;
    double[] assignmentProbs = null;
    int currentAssignment = 0;
    
    String[] baseString = {"A", "C", "G", "T"};
    
    Assembly[] assemblies = null; // set of possible data points
    int[][][] strandData = null; // [tp][strand][base] top two sets of reads on each strand
    int[][] strandCount = null; // [tp][strand] number of reads on each strand
    int[][] data = null; // [tp][base]
    int[] totData = null; // [tp]
    int sumTotData = 0; // sum of totData over all tp
    double[][] probStrand = null;
    double[][] logLikeAssignPoly = null;
    double[][] logStrandChoose = null; // N_s choose n_s for each [tp][strand]
    double[] logTotChoose = null;  // N choose n for each [tp];
    int[] activeBase = {-999, -999};  // active bases at that site
    int nActiveBase = 0;   // number of active bases
  
    boolean variable = false;
    boolean uncertain = false;
    boolean active = true;
    boolean[] missing = null;   // time point is missing
    
    boolean trial = false;
    
    double[][] alpha_loc = null;
    double alpha_e_loc = 0.0;
    double beta_e_loc = 0.0;
    double S_loc = 0.0;
    double F0_loc = 0.0;
    
    Site(int iSite, int nTimePoints, int nHaplo, int maxBases, int[] nAssignments, int[][] assign) {
        this.iSite = iSite;
        this.nTimePoints = nTimePoints;
        this.nHaplo = nHaplo;
        this.maxBases = maxBases;
        this.nAssignments = nAssignments;
        this.assign = assign;
        assemblies = new Assembly[nTimePoints]; 
        strandData = new int[nTimePoints][2][2];  
        data = new int[nTimePoints][2]; 
        totData = new int[nTimePoints];  
        strandCount = new int[nTimePoints][2];
        logStrandChoose = new double[nTimePoints][2]; 
        logTotChoose = new double[nTimePoints];
        missing = new boolean[nTimePoints];
        logLikeAssignPoly = new double[nTimePoints][nAssignments[maxBases]];
        probStrand = new double[nTimePoints][2];
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            missing[iTimePoint] = true; // initially assume data point doesn't exist
        }
    }
  
    void addAssembly(int iTimePoint, Assembly assembly) {  
        missing[iTimePoint] = false;
        assemblies[iTimePoint] = assembly;
    }
    
    boolean process() {
        identifyActive();
        if (!active) {
            return false;
        }
        if (uncertain) {
            return false;
        }
        makeSums();
        return ((Cluster.random.nextDouble() < Cluster.useFrac) && nActiveBase > 0);
    }
    
    void identifyActive() {  // identify which bases dominate
        double[] maxMinAmt = new double[4];    // accumulate minimal sums
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            if (!missing[iTimePoint] && assemblies[iTimePoint].hasData()) { 
                for (int iBase = 0; iBase < 4; iBase++) {
                    maxMinAmt[iBase] = Math.max(maxMinAmt[iBase], assemblies[iTimePoint].getMinAmt()[iBase]);
                }
            }
        }
        Vector<Double> maxMinAmtVector = new Vector<>();
        Hashtable<Double, Integer> maxMinAmtHash = new Hashtable<>();
        for (int iBase = 0; iBase < 4; iBase++) {
            if (maxMinAmt[iBase] > 1.0E-8) {
                maxMinAmtVector.add(maxMinAmt[iBase]);
                maxMinAmtHash.put(maxMinAmt[iBase], iBase);
            }
        }
        if (maxMinAmtVector.size() == 0) {
            active = false;
        } else {
            Collections.sort(maxMinAmtVector);
            nActiveBase = Math.min(maxBases, maxMinAmtVector.size());
            activeBase = new int[nActiveBase];
            for (int iPoint = 0; iPoint < nActiveBase; iPoint++) {
                int iBase = maxMinAmtVector.size() - 1 - iPoint;
                activeBase[iPoint] = maxMinAmtHash.get(maxMinAmtVector.get(iBase));
            } 
            
            if (maxMinAmtVector.size() > maxBases) {
                int iPoint = maxMinAmtVector.size() - 1 - maxBases;
                int aBase = maxMinAmtHash.get(maxMinAmtVector.get(iPoint));
                System.out.println(aBase);
                if (maxMinAmt[aBase] > 0.001) {
                    System.out.println("Site " + iSite + " Ignoring " 
                            + baseString[maxMinAmtHash.get(maxMinAmtVector.get(iPoint))] 
                            + "\t" + maxMinAmtVector.get(iPoint) + "\t" + maxMinAmtHash.get(maxMinAmtVector.get(iPoint)));
                    for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                        if (!missing[iTimePoint] && assemblies[iTimePoint].hasData()) { 
                            System.out.println(Arrays.toString(assemblies[iTimePoint].strandReads[0]) + " " + Arrays.toString(assemblies[iTimePoint].strandReads[1]));
                        }   
                    }
                    System.out.println();
                }
            }
        }
    }
    
    void makeSums() {
        for (int iTimePoint = 0; iTimePoint< nTimePoints; iTimePoint++) {
            if (!missing[iTimePoint] && assemblies[iTimePoint].hasData()) {
                for (int iStrand = 0; iStrand < 2; iStrand++) {
                    for (int iActive = 0; iActive < nActiveBase; iActive++) {
                        if (activeBase[iActive] >= 0) {
                            strandData[iTimePoint][iStrand][iActive] = 
                                assemblies[iTimePoint].strandReads[iStrand][activeBase[iActive]];
                            data[iTimePoint][iActive] += strandData[iTimePoint][iStrand][iActive];
                            totData[iTimePoint] += strandData[iTimePoint][iStrand][iActive];
                            sumTotData += strandData[iTimePoint][iStrand][iActive];
                        }
                    }
                }
                for (int iStrand = 0; iStrand < 2; iStrand++) {
                    logStrandChoose[iTimePoint][iStrand] = 
                            CombinatoricsUtils.factorialLog(strandData[iTimePoint][iStrand][0]+strandData[iTimePoint][iStrand][1])
                            - CombinatoricsUtils.factorialLog(strandData[iTimePoint][iStrand][0])
                            - CombinatoricsUtils.factorialLog(strandData[iTimePoint][iStrand][1]);
                }
                logTotChoose[iTimePoint] = CombinatoricsUtils.factorialLog(totData[iTimePoint])
                            - CombinatoricsUtils.factorialLog(data[iTimePoint][0])
                            - CombinatoricsUtils.factorialLog(data[iTimePoint][1]);
                if (data[iTimePoint][0] > 0) {
                    nActiveBase = Math.max(nActiveBase, 1);
                }
                if (data[iTimePoint][1] > 0) {
                    nActiveBase = Math.max(nActiveBase, 2);
                }

            }
        }
    }
    
    void setParams(int iStage, double alpha_e, double beta_e, double S, double F0, double[][] alpha){
        alpha_e_loc = alpha_e;
        beta_e_loc = beta_e;
        S_loc = S;
        F0_loc = F0;
        alpha_loc = alpha;
        if (iStage == 0) {
//            double[][] logLikeAssignPoly = new double[nTimePoints][nAssignments];
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                if (!missing[iTimePoint]){
                    for (int iAssignment = 1; iAssignment < nAssignments[maxBases]-1; iAssignment++) {
                        double[] sumAlpha = new double[2];
                        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                            sumAlpha[assign[iAssignment][iHaplo]] += alpha_loc[iTimePoint][iHaplo];
                        }
                        logLikeAssignPoly[iTimePoint][iAssignment] = 
                            logStrandChoose[iTimePoint][0] + logStrandChoose[iTimePoint][1]
                            + Beta.logBeta(data[iTimePoint][1]+sumAlpha[0], data[iTimePoint][0]+sumAlpha[1])
                            - Beta.logBeta(sumAlpha[0], sumAlpha[1]);
                    }
                }
            }
        }
    }
    
    
    double computeLogLikelihoodAlpha(double[] alpha, int iTimePoint) { // alpha[][] = iTP, iHaplo, iTP
        if (currentAssignment == 0 || currentAssignment == nAssignments[maxBases]-1 || missing[iTimePoint]) {
            return 0.0;
        }
        double totalProb = 0.0;
        double logProb = 0.0;           
        double[] sumAlpha = new double[2];
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
            sumAlpha[assign[currentAssignment][iHaplo]] += alpha[iHaplo];
        }
        logProb = Math.log(probStrand[iTimePoint][0]) + Math.log(probStrand[iTimePoint][1])
            + logStrandChoose[iTimePoint][0] + logStrandChoose[iTimePoint][1]
            + Beta.logBeta(data[iTimePoint][1]+sumAlpha[0], data[iTimePoint][0]+sumAlpha[1])
            - Beta.logBeta(sumAlpha[0], sumAlpha[1]);
        return logProb;
    }

        
    double computeLogLikelihood(double alpha_e, double beta_e, double S, double F0, boolean printHaplo) { // alpha[][] = iHaplo, iTP
        double[] logLikeAssign = new double[nAssignments[maxBases]];  // log likelihood of the site for each assignment
        double totalProb = 0.0;
        assignmentProbs = new double[nAssignments[maxBases]];
        double[][] likeAssignSNoS = new double[nTimePoints][2];
        for (int iAssignment = 0; iAssignment < nAssignments[maxBases]; iAssignment++) {
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                if (missing[iTimePoint]) {
                    likeAssignSNoS[iTimePoint][0] = 1.0;
                    likeAssignSNoS[iTimePoint][1] = 0.0;
                }
                if (!missing[iTimePoint]) {
                    double[][][] probStrandSNoS = new double[nAssignments[maxBases]][2][2];
                    
                    if (iAssignment == 0) {
                        for (int iStrand = 0; iStrand < 2; iStrand++) {
                            double logProbTerm = logStrandChoose[iTimePoint][iStrand]
                            + Beta.logBeta(strandData[iTimePoint][iStrand][1]+alpha_e, strandData[iTimePoint][iStrand][0]+beta_e)
                            - Beta.logBeta(alpha_e, beta_e);
                            probStrandSNoS[iAssignment][iStrand][0] = (1.0-S)*Math.exp(logProbTerm);
                            probStrandSNoS[iAssignment][iStrand][1] = S/(strandData[iTimePoint][iStrand][0]+strandData[iTimePoint][iStrand][1]+1.0);
                            probStrand[iTimePoint][iStrand] = probStrandSNoS[iAssignment][iStrand][0] 
                                    + probStrandSNoS[iAssignment][iStrand][1];
                        }
                        logLikeAssign[iAssignment] += Math.log(probStrand[iTimePoint][0]) + Math.log(probStrand[iTimePoint][1]);
                        likeAssignSNoS[iTimePoint][0] += probStrandSNoS[iAssignment][0][0]*probStrandSNoS[iAssignment][1][0];
                        likeAssignSNoS[iTimePoint][1] += probStrand[iTimePoint][0]*probStrand[iTimePoint][1]
                                - probStrandSNoS[iAssignment][0][0]*probStrandSNoS[iAssignment][1][0];
                    } else if (iAssignment == nAssignments[maxBases] - 1) {
                        for (int iStrand = 0; iStrand < 2; iStrand++) {
                            double logProbTerm = logStrandChoose[iTimePoint][iStrand]
                                + Beta.logBeta(strandData[iTimePoint][iStrand][0]+alpha_e, strandData[iTimePoint][iStrand][1]+beta_e)
                                - Beta.logBeta(alpha_e, beta_e);
                            probStrandSNoS[iAssignment][iStrand][0] = (1.0 - S) * Math.exp(logProbTerm);
                            probStrandSNoS[iAssignment][iStrand][1] = S/(strandData[iTimePoint][iStrand][0] + strandData[iTimePoint][iStrand][1]+1.0);
                            probStrand[iTimePoint][iStrand] = probStrandSNoS[iAssignment][iStrand][0] 
                                    + probStrandSNoS[iAssignment][iStrand][1];
                        }
                        logLikeAssign[iAssignment] += Math.log(probStrand[iTimePoint][0]) + Math.log(probStrand[iTimePoint][1]);
                        likeAssignSNoS[iTimePoint][0] += probStrandSNoS[iAssignment][0][0]*probStrandSNoS[iAssignment][1][0];
                        likeAssignSNoS[iTimePoint][1] += probStrand[iTimePoint][0]*probStrand[iTimePoint][1]
                                - probStrandSNoS[iAssignment][0][0]*probStrandSNoS[iAssignment][1][0];

                    } else {
                        for (int iStrand = 0; iStrand < 2; iStrand++) {
                            double logProbTerm = 
                                Beta.logBeta(alpha_e, strandData[iTimePoint][iStrand][0]+strandData[iTimePoint][iStrand][1] + beta_e) 
                                    - Beta.logBeta(alpha_e, beta_e);
                             probStrand[iTimePoint][iStrand] = (1.0-S)*Math.exp(logProbTerm)
                                     + S/(strandData[iTimePoint][iStrand][0]+strandData[iTimePoint][iStrand][1]+1.0);
                        }
                        logLikeAssign[iAssignment] += Math.log(probStrand[iTimePoint][0]) + Math.log(probStrand[iTimePoint][1])
                                + logLikeAssignPoly[iTimePoint][iAssignment];
                        likeAssignSNoS[iTimePoint][0] += probStrand[iTimePoint][0]*probStrand[iTimePoint][1]*Math.exp(logLikeAssignPoly[iTimePoint][iAssignment]);
                    }   // end of elses
                }  // end of conditional on timepoint existing
            }   // end of loop over timepoints
 
            if (iAssignment == 0 || iAssignment == nAssignments[maxBases] - 1) {
                assignmentProbs[iAssignment] = (1.0 - F0) * 0.5 * Math.exp(logLikeAssign[iAssignment]);
            } else {
                assignmentProbs[iAssignment] += (F0 / (nAssignments[maxBases] - 2.0)) * Math.exp(logLikeAssign[iAssignment]);
            }
            totalProb += assignmentProbs[iAssignment];
        } // end of loop over assignments
        currentAssignment = iSelect(assignmentProbs);
        if (printHaplo) {
            double maxNoSProb = 1.0;
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                maxNoSProb *= likeAssignSNoS[iTimePoint][0]/(likeAssignSNoS[iTimePoint][0]+likeAssignSNoS[iTimePoint][1]);
            }
            System.out.print(iSite);
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                double[] sumProb = new double[2];
                for (int iAssignment = 0; iAssignment < nAssignments[maxBases]; iAssignment++) {
                    sumProb[assign[iAssignment][iHaplo]] += assignmentProbs[iAssignment];
                }
                if (sumProb[0] > sumProb[1]) {
                    System.out.format("\t%s\t%5.3f", baseString[activeBase[0]], maxNoSProb*sumProb[0]);
                } else {
                    System.out.format("\t%s\t%5.3f", baseString[activeBase[1]], maxNoSProb*sumProb[1]);
                }
            }
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                if (iTimePoint == 0) {
                    System.out.print("\t");
                } else {
                    System.out.print(" ");
                }
                System.out.print(Arrays.toString(data[iTimePoint]));
            }
            System.out.println();
        }
        return Math.log(totalProb);
    }
    
    int iSelect(double[] probs) {   // Selects integer from a specified array of probabilities
        int nPoss = probs.length;
        double summ = 0.0;
        for (int iPoss = 0; iPoss < nPoss; iPoss++) {
            summ += probs[iPoss];
        }
        for (int iPoss = 0; iPoss < nPoss; iPoss++) {
            probs[iPoss] /= summ;
        }
        double randomNumber = Cluster.random.nextDouble();
        for (int iPoss = 0; iPoss < nPoss-1; iPoss++) {
            if (probs[iPoss] > randomNumber) {
                return iPoss;
            } else {
                randomNumber -= probs[iPoss];
            }
        }
        return (nPoss-1);
    }

}
