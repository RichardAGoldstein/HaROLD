/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Cluster;

import java.util.Arrays;
import org.apache.commons.math3.special.Gamma;

/**
 *
 * @author rgoldst
 */
public class Assignment {
    
    int[] assign = null;
    boolean[] present = new boolean[4];
    int nHaplo = 0;
    int nPresent = 0;
    int nAbsent = 0;
    int nTimePoints = 0;
    
    double[][] alphaHep = null;
    double[] quartiles = null;
    int nQuartiles = 0;
    double S = 0.0;
    double[][] alpha_nuc = null;
    double alpha_C = 0.0;
    double alpha_E = 0.0;
    double bEAlpha_E = 0.0;
    double[] sumAlpha = null;
    double[][] logGammaAlpha_nuc = null;
    double[] logGammaSumAlpha = null;
    double logGammaAlpha_E = 0.0;
    double logGammaAlpha_C = 0.0;
    double logGammaSumCPlusE = 0.0;
    
    Assignment(int iAssign, int nHaplo) {
        this.nHaplo = nHaplo;
        assign = new int[nHaplo];
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {    // Loop over possible haplotypes
            assign[iHaplo] = (iAssign / pow(Cluster.maxBases, iHaplo)) % (Cluster.maxBases);
            present[assign[iHaplo]] = true;
        }       
        for (int iBase = 0; iBase < 4; iBase++) {
            if (present[iBase]) {
                nPresent++;
            }
        }
        nAbsent = 4 - nPresent;
    }
    
    void setAlphas(double[][] alphaHap, double[] quartiles, double S) {
        this.alphaHep = alphaHep;
        this.quartiles = quartiles;
        nQuartiles = quartiles.length;
        this.S = S;
        
        nTimePoints = alphaHap.length;
        alpha_nuc = new double[nTimePoints][4];       
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                alpha_nuc[iTimePoint][assign[iHaplo]] += alphaHap[iTimePoint][iHaplo];
            }
        }
    }

    
    double computeAssignmentLogLikelihood(int iTimePoint, int[][] strandReads, int[] reads, int[] totStrand, boolean siteConserved ) {
        double[] logFitness1 = new double[2];
        double[] logFitness2 = new double[2];
        logFitness2[0] = Gamma.logGamma(4.0) - Gamma.logGamma(totStrand[0]+4.0);
        logFitness2[1] = Gamma.logGamma(4.0) - Gamma.logGamma(totStrand[1]+4.0);
        for (int iBase = 0; iBase < 4; iBase++) {
            if (reads[iBase] > 0) {
                double[][] terms = new double[2][nQuartiles];
                double maxTerm = -1.0E10;
                for (int iQuartile = 0; iQuartile < nQuartiles; iQuartile++) {
                    double prob = alpha_nuc[iTimePoint][iBase] + (1.0 - 4.0 * alpha_nuc[iTimePoint][iBase]) * quartiles[iQuartile];
                    terms[0][iQuartile] = strandReads[0][iBase] * Math.log(prob) - Math.log(1.0*nQuartiles);
                    terms[1][iQuartile] = strandReads[1][iBase] * Math.log(prob) - Math.log(1.0*nQuartiles);
                    maxTerm = Math.max(Math.max(maxTerm, terms[0][iQuartile]), terms[1][iQuartile]);
                }
                double[] sumTerms = new double[2];
                for (int iQuartile = 0; iQuartile < nQuartiles; iQuartile++) {
                    sumTerms[0] += Math.exp(terms[0][iQuartile] - maxTerm);
                    sumTerms[1] += Math.exp(terms[1][iQuartile] - maxTerm);
                }    
                logFitness1[0] += maxTerm + Math.log(sumTerms[0]);
                logFitness1[1] += maxTerm + Math.log(sumTerms[1]);
                logFitness2[0] += Gamma.logGamma(strandReads[0][iBase] + 1.0);
                logFitness2[1] += Gamma.logGamma(strandReads[1][iBase] + 1.0);
            }
        }
        double logFit11 = logFitness1[0] + logFitness1[1];
        double logFit12 = logFitness1[0] + logFitness2[1];
        double logFit21 = logFitness2[0] + logFitness1[1];
        double logFit22 = logFitness2[0] + logFitness2[1];
        double logFitnessMax = Math.max(Math.max(logFit11, logFit12), Math.max(logFit21, logFit22));
        double logFitness = (1.0 - S) * (1.0 - S) * Math.exp(logFit11 - logFitnessMax)
                            + (1.0 - S) * S * Math.exp(logFit12 - logFitnessMax)
                            + (1.0 - S) * S * Math.exp(logFit21 - logFitnessMax)
                            + S * S * Math.exp(logFit22 - logFitnessMax);
        logFitness = logFitnessMax + Math.log(logFitness);
//        System.out.println(iTimePoint + "\t" + Arrays.toString(strandReads[0]) + "  " 
//                + Arrays.toString(strandReads[1]) + "\t" + 
//                Arrays.toString(assign) + "\t" + Arrays.toString(alpha_nuc[iTimePoint]) + "\t" + logFitness
//                + "\t" + logFit11 + "  " + logFit12 + "  " + logFit21 + "  " + logFit22);
        return logFitness;
    }
    
    
    
    int pow (int a, int b) {  // Computes powers
        if ( b == 0)     return 1;
        if ( b == 1)     return a;
        if (b%2 == 0)    return     pow ( a * a, b/2); //even a=(a^2)^b/2
        else             return a * pow ( a * a, (b-1)/2); //odd  a=a*(a^2)^b/2
    }

    
}
