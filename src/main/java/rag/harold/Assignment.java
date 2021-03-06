package rag.harold;

import java.util.Arrays;

/**
 * @author rgoldst
 */
public class Assignment {

    private final GammaCalc gamma;
    private final boolean verbose;
    int[] assign = null;
    boolean[] presentBase = new boolean[4];
    int nPresent = 0;
    private int nHaplo = 0;
    private int nAbsent = 0;
    private int nTimePoints = 0;
    private double[][] currentPiNuc = null;
    private double[][] currentAlphaObs = null;
    private double[] currentSumAlphaObs = null;
    private double currentAlpha0 = 0.0;
    private double currentAlphaE = 0.0;

    Assignment(int iAssign, int nHaplo, GammaCalc gammaCalc, boolean verbose) {
        this.gamma = gammaCalc;
        this.nHaplo = nHaplo;
        this.verbose = verbose;
        assign = new int[nHaplo];
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {    // Loop over possible haplotypes
            assign[iHaplo] = (iAssign / pow(Constants.MAX_BASES, iHaplo)) % (Constants.MAX_BASES);
            presentBase[assign[iHaplo]] = true;
        }
        for (int iBase = 0; iBase < 4; iBase++) {
            if (presentBase[iBase]) {
                nPresent++;
            }
        }
        nAbsent = 4 - nPresent;
        if (this.verbose) {
            System.out.println(iAssign + "\t" + Arrays.toString(assign));
        }
    }

    void setAllParams(double[][] piHap, double[] alphaParams) {
        nTimePoints = piHap.length;
        currentPiNuc = new double[nTimePoints][4];
        currentAlphaObs = new double[nTimePoints][4];
        currentSumAlphaObs = new double[nTimePoints];
        currentAlpha0 = alphaParams[0] * (1.0 - alphaParams[1]) / alphaParams[1];
        currentAlphaE = (1.0 - alphaParams[0]) * (1.0 - alphaParams[1]) / alphaParams[1];
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            setPiHap(iTimePoint, piHap[iTimePoint]);
        }
    }

    void setSinglePiHap(int iTimePoint, double[] piHap) {
        Arrays.fill(currentPiNuc[iTimePoint], 0.0);
        Arrays.fill(currentAlphaObs[iTimePoint], 0.0);
        currentSumAlphaObs[iTimePoint] = 0.0;
        setPiHap(iTimePoint, piHap);
    }

    void setPiHap(int iTimePoint, double[] piHap) {
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
            currentPiNuc[iTimePoint][assign[iHaplo]] += piHap[iHaplo];
        }
        for (int iBase = 0; iBase < 4; iBase++) {
            currentAlphaObs[iTimePoint][iBase]
                    = currentPiNuc[iTimePoint][iBase] * currentAlpha0
                    + (1.0 - currentPiNuc[iTimePoint][iBase]) * currentAlphaE;
            currentSumAlphaObs[iTimePoint] += currentAlphaObs[iTimePoint][iBase];
        }
    }

    double computeAssignmentLogLikelihood(int iTimePoint, int[][] strandReads, int[] reads, int[] totStrand, boolean siteConserved) {
        double[] logLikelihoodStrand = new double[2];
        double g1 = this.gamma.logGamma(currentSumAlphaObs[iTimePoint]);
        double[] currentAlphaObs_iTimePoint = currentAlphaObs[iTimePoint];

        for (int iStrand = 0; iStrand < 2; iStrand++) {
            logLikelihoodStrand[iStrand] = g1 - this.gamma.logGamma(currentSumAlphaObs[iTimePoint] + totStrand[iStrand]);
            for (int iBase = 0; iBase < 4; iBase++) {
                if (strandReads[iStrand][iBase] > 0) {
                    logLikelihoodStrand[iStrand] += this.gamma.logGamma(currentAlphaObs_iTimePoint[iBase] + strandReads[iStrand][iBase])
                            - this.gamma.logGamma(currentAlphaObs_iTimePoint[iBase]);
                }
            }
        }
        double logFitness = logLikelihoodStrand[0] + logLikelihoodStrand[1];
        return logFitness;
    }


    int pow(int a, int b) {  // Computes powers
        if (b == 0) return 1;
        if (b == 1) return a;
        if (b % 2 == 0) return pow(a * a, b / 2); //even a=(a^2)^b/2
        else return a * pow(a * a, (b - 1) / 2); //odd  a=a*(a^2)^b/2
    }


}
