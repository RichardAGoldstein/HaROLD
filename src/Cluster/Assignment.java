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
    
    void setAlphas(double[][] alphaHap, double alpha_C, double alpha_E) {
        this.alpha_C = alpha_C;
        this.alpha_E = alpha_E;
//        this.bEAlpha_E = nAbsent * alpha_E;
        this.bEAlpha_E = alpha_E;         // Removing b term
        
        nTimePoints = alphaHap.length;
        alpha_nuc = new double[nTimePoints][4];
        logGammaAlpha_nuc = new double[nTimePoints][4];
        sumAlpha = new double[nTimePoints];
        logGammaSumAlpha = new double[nTimePoints];
        
        logGammaAlpha_E = Gamma.logGamma(alpha_E);
        logGammaAlpha_C = Gamma.logGamma(alpha_C);
        logGammaSumCPlusE = Gamma.logGamma(alpha_C + bEAlpha_E);
        
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                alpha_nuc[iTimePoint][assign[iHaplo]] += alphaHap[iTimePoint][iHaplo];
            }

            for (int iNuc = 0; iNuc < 4; iNuc++) {
                sumAlpha[iTimePoint] += alpha_nuc[iTimePoint][iNuc];
                logGammaAlpha_nuc[iTimePoint][iNuc] = Gamma.logGamma(alpha_nuc[iTimePoint][iNuc]);
            }
            logGammaSumAlpha[iTimePoint] = Gamma.logGamma(sumAlpha[iTimePoint]);
        }

    }

    
    double computeAssignmentLogLikelihood(int iTimePoint, int[][] strandReads, int[] reads, int[] totStrand, boolean siteConserved ) {
        int[] N_Cs = new int[2];
        int N_C = 0;
        int[] errors = new int[2];
        double logFitness = 2.0 * (logGammaSumCPlusE - logGammaAlpha_C);

        for (int iBase = 0; iBase < 4; iBase++) {
            if (present[iBase]) {
                if (reads[iBase] > 0) {
                    logFitness += Gamma.logGamma(reads[iBase] + alpha_nuc[iTimePoint][iBase]) - logGammaAlpha_nuc[iTimePoint][iBase];
                    N_Cs[0] += strandReads[0][iBase];
                    N_Cs[1] += strandReads[1][iBase];
                    N_C += reads[iBase];
                }
            } else {
                errors[0] += strandReads[0][iBase];
                errors[1] += strandReads[1][iBase];
            }
        }
        logFitness += (logGammaSumAlpha[iTimePoint] - Gamma.logGamma(N_C + sumAlpha[iTimePoint]))
                + (Gamma.logGamma(N_Cs[0] + alpha_C) - Gamma.logGamma(totStrand[0] + alpha_C  + bEAlpha_E))
                + (Gamma.logGamma(N_Cs[1] + alpha_C) - Gamma.logGamma(totStrand[1] + alpha_C  + bEAlpha_E));
        if (errors[0] > 0) {
            logFitness += Gamma.logGamma(errors[0] + alpha_E) - logGammaAlpha_E;
        }
        if (errors[1] > 0) {
            logFitness += Gamma.logGamma(errors[1] + alpha_E) - logGammaAlpha_E;
        }
        if (false) {
            System.out.println(logFitness + "\t" + N_C + "\t" + N_Cs[0] + "\t" 
                    + N_Cs[1] + "\t" + Arrays.toString(reads) + "\t" + Arrays.toString(assign) + "\t" 
                    + Arrays.toString(alpha_nuc));
        }
        return logFitness;
    }
    
    
    
    int pow (int a, int b) {  // Computes powers
        if ( b == 0)     return 1;
        if ( b == 1)     return a;
        if (b%2 == 0)    return     pow ( a * a, b/2); //even a=(a^2)^b/2
        else             return a * pow ( a * a, (b-1)/2); //odd  a=a*(a^2)^b/2
    }

    
}
