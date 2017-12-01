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
    boolean[] presentBase = new boolean[4];
    int nHaplo = 0;
    int nPresent = 0;
    int nAbsent = 0;
    int nTimePoints = 0;
    
    double[][] piHap = null;
    double[][] piNuc = null;
    double[][] alphaObs = null;
    double alpha0 = 0.0;
    double alphaE = 0.0;
    
    Assignment(int iAssign, int nHaplo) {
        this.nHaplo = nHaplo;
        assign = new int[nHaplo];
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {    // Loop over possible haplotypes
            assign[iHaplo] = (iAssign / pow(Cluster.maxBases, iHaplo)) % (Cluster.maxBases);
            presentBase[assign[iHaplo]] = true;
        }       
        for (int iBase = 0; iBase < 4; iBase++) {
            if (presentBase[iBase]) {
                nPresent++;
            }
        }
        nAbsent = 4 - nPresent;
    }
    
    void setParams(double[][] piHap, double alpha0, double alphaE) {
        setPiHap(piHap);
        setAlphas(alpha0, alphaE);
    }

    void setParams(double[][] piHap, double[] alphaParams) {
        setPiHap(piHap);
        setAlphas(alphaParams);
    }

    
    void setAlphas(double alpha0, double alphaE) {
        this.alpha0 = alpha0;
        this.alphaE = alphaE;
    }

    void setAlphas(double[] alphaParams) {
        this.alpha0 = alphaParams[0];
        this.alphaE = alphaParams[1];
    }
    
    void setPiHap(double[][] piHap) {
        this.piHap = piHap;
        nTimePoints = piHap.length;
        piNuc = new double[nTimePoints][4];
        alphaObs = new double[nTimePoints][4];
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                piNuc[iTimePoint][assign[iHaplo]] += piHap[iTimePoint][iHaplo];
            }
            for (int iBase = 0; iBase < 4; iBase++) {
                alphaObs[iTimePoint][iBase] = piNuc[iTimePoint][iBase] * alpha0 - (1.0 - piNuc[iTimePoint][iBase]) * alphaE;
            }
        }
    }    
    
    
    void setPiHap(int iTimePoint, double[] piHap) {
        piNuc = new double[nTimePoints][4];
        alphaObs = new double[nTimePoints][4];
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
            piNuc[iTimePoint][assign[iHaplo]] += piHap[iHaplo];
        }
        for (int iBase = 0; iBase < 4; iBase++) {
            this.piHap[iTimePoint][iBase] = piHap[iBase];
            alphaObs[iTimePoint][iBase] = piNuc[iTimePoint][iBase] * alpha0 - (1.0 - piNuc[iTimePoint][iBase]) * alphaE;
        }
    } 
    
    double computeAssignmentLogLikelihood(int iTimePoint, int[][] strandReads, int[] reads, int[] totStrand, boolean siteConserved ) {
        double[] logLikelihoodStrand = new double[2];
        if (siteConserved) {
            for (int iStrand = 0; iStrand < 2; iStrand++) {
                logLikelihoodStrand[iStrand] = Gamma.logGamma(alpha0 + 3.0 * alphaE)
                        - Gamma.logGamma(alpha0 + 3.0 * alphaE + totStrand[iStrand])
                        + Gamma.logGamma(alpha0 + totStrand[iStrand])
                        - Gamma.logGamma(alpha0);
            }
        } else {
            for (int iStrand = 0; iStrand < 2; iStrand++) {
                logLikelihoodStrand[iStrand] = Gamma.logGamma(alpha0 + 3.0 * alphaE)
                        - Gamma.logGamma(alpha0 + 3.0 * alphaE + totStrand[iStrand]);
                for (int iBase = 0; iBase < 4; iBase++) {
                    if (strandReads[iStrand][iBase] > 0) {
                        logLikelihoodStrand[iStrand] += Gamma.logGamma(alphaObs[iTimePoint][iBase] + strandReads[iStrand][iBase])
                                - Gamma.logGamma(alphaObs[iTimePoint][iBase]);
                    }
                }
            }          
        }
        double logFitness = logLikelihoodStrand[0] + logLikelihoodStrand[1];
        return logFitness;
    }
    
    
    
    int pow (int a, int b) {  // Computes powers
        if ( b == 0)     return 1;
        if ( b == 1)     return a;
        if (b%2 == 0)    return     pow ( a * a, b/2); //even a=(a^2)^b/2
        else             return a * pow ( a * a, (b-1)/2); //odd  a=a*(a^2)^b/2
    }

    
}
