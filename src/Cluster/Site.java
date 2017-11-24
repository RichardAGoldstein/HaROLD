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
 * Container for the reads at a single site
 * 
 * @author rgoldst
 */
public class Site {   
    int iSite;
    int nTimePoints;
    int nHaplo = 0;
    int nBases = 0;
    int[][] assign = null;
    int nAssignments = 0;
    
    
    double[] assignmentProbs = null;
    int currentAssignment = 0;
    
    String[] baseString = {"A", "C", "G", "T"};
    
    
    int[][][] strandReads = null; // [tp][strand][base] top two sets of reads on each strand
    int[][] totStrand = null; // [tp][strand] number of reads on each strand
    int[][] reads = null; // [tp][base]
    int[] totReads = null; // [tp]
    int[] timePointConservedBase = null; 
    
    
    
    double[][] probStrand = null;
    double[][] logLikeAssignPoly = null;
    double[][] logStrandChoose = null; // N_s choose n_s for each [tp][strand]
    double[] logTotChoose = null;  // N choose n for each [tp];
    
    boolean siteActive = false;
    int nPresentBase = 0;   // number of present bases
    boolean[] presentBase = new boolean[4];
    boolean siteConserved = false;
    boolean[] timePointConserved = null;
    boolean[] timePointHasData = null;

    
    double[][] alpha_loc = null;
    double alpha_e_loc = 0.0;
    double beta_e_loc = 0.0;
    double S_loc = 0.0;
    double F0_loc = 0.0;
    
    Site(int iSite, int nTimePoints, int nHaplo, int[][] assign) {
        this.iSite = iSite;
        this.nTimePoints = nTimePoints;
        this.nHaplo = nHaplo;
        this.assign = assign;
        nAssignments = assign.length;
        strandReads = new int[nTimePoints][2][4]; // [tp][strand][base] top two sets of reads on each strand
        totStrand = new int[nTimePoints][2]; // [tp][strand] number of reads on each strand
        reads = new int[nTimePoints][4]; // [tp][base]
        totReads = new int[nTimePoints]; // [tp]
        timePointConserved = new boolean[nTimePoints];
        timePointHasData = new boolean[nTimePoints];
        timePointConservedBase = new int[nTimePoints];
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            timePointHasData[iTimePoint] = false; // assume timepoint does not have data
            timePointConserved[iTimePoint] = false; // assume timepoint is not conserved
        }
    }
  
    void addTimePoint(int iTimePoint, String line) { 
        String[] words = line.split(",");
        for (int iBase = 0; iBase < 4; iBase++) {   // compute various sums of reads
            for (int iStrand = 0; iStrand < 2; iStrand++) {
                strandReads[iTimePoint][iStrand][iBase]=Integer.parseInt(words[(2*iBase)+iStrand+3]);
                totStrand[iTimePoint][iStrand] += strandReads[iTimePoint][iStrand][iBase];
                reads[iTimePoint][iBase] += strandReads[iTimePoint][iStrand][iBase];
                totReads[iTimePoint] += strandReads[iTimePoint][iStrand][iBase];
            }
            if (reads[iTimePoint][iBase] > 0) {
                presentBase[iBase] = true;
            }
        }
        if (totReads[iTimePoint] > 0) {   // has data
            timePointHasData[iTimePoint] = true;
        }
        if (Math.max(Math.max(reads[iTimePoint][0],reads[iTimePoint][1]),
                Math.max(reads[iTimePoint][2],reads[iTimePoint][3])) == totReads[iTimePoint]) {  // timePointConserved site
            timePointConserved[iTimePoint] = true;
            for (int iBase = 0; iBase < 4; iBase++) {
                if (reads[iTimePoint][iBase] == totReads[iTimePoint]) {
                    timePointConservedBase[iTimePoint] = iBase;
                }
            }
        }
    }
    
    double computeLogFitness(int iTimePoint, double[] alpha_nuc, double alpha_C, double alpha_E) {
        double logFitness = 0.0;
        double bAlpha_E = 0.0;
        double alphaSum = 0.0;
        double[] N_Cs = new double[2];
        for (int iBase = 0; iBase < 4; iBase++) {
            if (alpha_nuc[iBase] < 1.0E-10) {
                bAlpha_E += alpha_E;
                logFitness += Gamma.logGamma(strandReads[iTimePoint][0][iBase] + alpha_E)
                        + Gamma.logGamma(strandReads[iTimePoint][1][iBase] + alpha_E); 
            } else {
                alphaSum += alpha_nuc[iBase];
                N_Cs[0] += strandReads[iTimePoint][0][iBase];
                N_Cs[1] += strandReads[iTimePoint][1][iBase];
                logFitness += Gamma.logGamma(reads[iTimePoint][iBase] + alpha_nuc[iBase]);
            }
        }
        double N_C = N_Cs[0] + N_Cs[1];
        logFitness += Gamma.logGamma(N_Cs[0] + alpha_C) + Gamma.logGamma(N_Cs[1] + alpha_C) - 
                (Gamma.logGamma(N_C + alphaSum) + Gamma.logGamma(totStrand[iTimePoint][0] + alpha_C + bAlpha_E) 
                    + Gamma.logGamma(totStrand[iTimePoint][1] + alpha_C + bAlpha_E));
        return logFitness;
    }
    
    
    
    boolean isActive() {
        siteActive = false;
        for (int iBase = 0; iBase < 4; iBase++) {
            if (presentBase[iBase]) {
                nPresentBase++;
                siteActive = true;
            }
        }
        return siteActive;
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
                if (timePointHasData[iTimePoint]){
                    for (int iAssignment = 1; iAssignment < nAssignments-1; iAssignment++) {
                        double[] sumAlpha = new double[2];
                        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                            sumAlpha[assign[iAssignment][iHaplo]] += alpha_loc[iTimePoint][iHaplo];
                        }
                        logLikeAssignPoly[iTimePoint][iAssignment] = 
                            logStrandChoose[iTimePoint][0] + logStrandChoose[iTimePoint][1]
                            + Beta.logBeta(reads[iTimePoint][1]+sumAlpha[0], reads[iTimePoint][0]+sumAlpha[1])
                            - Beta.logBeta(sumAlpha[0], sumAlpha[1]);
                    }
                }
            }
        }
    }
    
    
    double computeLogLikelihoodAlpha(double[] alpha, int iTimePoint) { // alpha[][] = iTP, iHaplo, iTP
        if (currentAssignment == 0 || currentAssignment == nAssignments-1 || !timePointHasData[iTimePoint]) {
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
            + Beta.logBeta(reads[iTimePoint][1]+sumAlpha[0], reads[iTimePoint][0]+sumAlpha[1])
            - Beta.logBeta(sumAlpha[0], sumAlpha[1]);
        return logProb;
    }

        
    double computeLogLikelihood(double alpha_e, double beta_e, double S, double F0, boolean printHaplo) { // alpha[][] = iHaplo, iTP
        double[] logLikeAssign = new double[nAssignments];  // log likelihood of the site for each assignment
        double totalProb = 0.0;
        assignmentProbs = new double[nAssignments];
        double[][] likeAssignSNoS = new double[nTimePoints][2];
        for (int iAssignment = 0; iAssignment < nAssignments; iAssignment++) {
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                if (!timePointHasData[iTimePoint]) {
                    likeAssignSNoS[iTimePoint][0] = 1.0;
                    likeAssignSNoS[iTimePoint][1] = 0.0;
                }
                if (timePointHasData[iTimePoint]) {
                    double[][][] probStrandSNoS = new double[nAssignments][2][2];
                    
                    if (iAssignment == 0) {
                        for (int iStrand = 0; iStrand < 2; iStrand++) {
                            double logProbTerm = logStrandChoose[iTimePoint][iStrand]
                            + Beta.logBeta(strandReads[iTimePoint][iStrand][1]+alpha_e, strandReads[iTimePoint][iStrand][0]+beta_e)
                            - Beta.logBeta(alpha_e, beta_e);
                            probStrandSNoS[iAssignment][iStrand][0] = (1.0-S)*Math.exp(logProbTerm);
                            probStrandSNoS[iAssignment][iStrand][1] = S/(strandReads[iTimePoint][iStrand][0]+strandReads[iTimePoint][iStrand][1]+1.0);
                            probStrand[iTimePoint][iStrand] = probStrandSNoS[iAssignment][iStrand][0] 
                                    + probStrandSNoS[iAssignment][iStrand][1];
                        }
                        logLikeAssign[iAssignment] += Math.log(probStrand[iTimePoint][0]) + Math.log(probStrand[iTimePoint][1]);
                        likeAssignSNoS[iTimePoint][0] += probStrandSNoS[iAssignment][0][0]*probStrandSNoS[iAssignment][1][0];
                        likeAssignSNoS[iTimePoint][1] += probStrand[iTimePoint][0]*probStrand[iTimePoint][1]
                                - probStrandSNoS[iAssignment][0][0]*probStrandSNoS[iAssignment][1][0];
                    } else if (iAssignment == nAssignments - 1) {
                        for (int iStrand = 0; iStrand < 2; iStrand++) {
                            double logProbTerm = logStrandChoose[iTimePoint][iStrand]
                                + Beta.logBeta(strandReads[iTimePoint][iStrand][0]+alpha_e, strandReads[iTimePoint][iStrand][1]+beta_e)
                                - Beta.logBeta(alpha_e, beta_e);
                            probStrandSNoS[iAssignment][iStrand][0] = (1.0 - S) * Math.exp(logProbTerm);
                            probStrandSNoS[iAssignment][iStrand][1] = S/(strandReads[iTimePoint][iStrand][0] + strandReads[iTimePoint][iStrand][1]+1.0);
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
                                Beta.logBeta(alpha_e, strandReads[iTimePoint][iStrand][0]+strandReads[iTimePoint][iStrand][1] + beta_e) 
                                    - Beta.logBeta(alpha_e, beta_e);
                             probStrand[iTimePoint][iStrand] = (1.0-S)*Math.exp(logProbTerm)
                                     + S/(strandReads[iTimePoint][iStrand][0]+strandReads[iTimePoint][iStrand][1]+1.0);
                        }
                        logLikeAssign[iAssignment] += Math.log(probStrand[iTimePoint][0]) + Math.log(probStrand[iTimePoint][1])
                                + logLikeAssignPoly[iTimePoint][iAssignment];
                        likeAssignSNoS[iTimePoint][0] += probStrand[iTimePoint][0]*probStrand[iTimePoint][1]*Math.exp(logLikeAssignPoly[iTimePoint][iAssignment]);
                    }   // end of elses
                }  // end of conditional on timepoint existing
            }   // end of loop over timepoints
 
            if (iAssignment == 0 || iAssignment == nAssignments - 1) {
                assignmentProbs[iAssignment] = (1.0 - F0) * 0.5 * Math.exp(logLikeAssign[iAssignment]);
            } else {
                assignmentProbs[iAssignment] += (F0 / (nAssignments - 2.0)) * Math.exp(logLikeAssign[iAssignment]);
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
                for (int iAssignment = 0; iAssignment < nAssignments; iAssignment++) {
                    sumProb[assign[iAssignment][iHaplo]] += assignmentProbs[iAssignment];
                }
            }
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                if (iTimePoint == 0) {
                    System.out.print("\t");
                } else {
                    System.out.print(" ");
                }
                System.out.print(Arrays.toString(reads[iTimePoint]));
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
