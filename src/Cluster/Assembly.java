/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Cluster;

import java.util.*;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.*;

/**
 * Container for data for single data point
 * Precomputes stuff for calculation of probPoly
 * 
 * @author rgoldst
 */
public class Assembly {
    String inputLine;
    int iSite = 0;
    int[][] strandReads = new int[2][4];  // strand and base
    int[] reads = new int[4];  // sum for each base
    int totReads = 0;
    int[] totStrand = new int[2];
    double[] minAmt = new double[4];
    
    boolean hasData = true;
    boolean conserved = false;
    
    Assembly(String inputLine){   // create data point and read in line of input file
        this.inputLine = inputLine;
        String[] words = inputLine.split(",");
        iSite = Integer.parseInt(words[0]);

        for (int iBase = 0; iBase < 4; iBase++) {
            for (int iStrand = 0; iStrand < 2; iStrand++) {
                strandReads[iStrand][iBase]=Integer.parseInt(words[(2*iBase)+iStrand+3]);
                totStrand[iStrand] += strandReads[iStrand][iBase];
                reads[iBase] += strandReads[iStrand][iBase];
                totReads += strandReads[iStrand][iBase];
            }
        }
        if (totReads == 0) {
            hasData = false;

        } else if (Math.max(Math.max(reads[0],reads[1]),Math.max(reads[2],reads[3])) == totReads) {
            conserved = true;
        } else {
            double topRatio = 0.00000;
            int iConsensus = -999;
            boolean printThis = false;
            for (int iBase = 0; iBase < 4; iBase++) {
                if (reads[iBase] > 0) {
                    int[] n = new int[2];
                    int[] m = new int[2];
                    double[] ratio = new double[2];
                    n[0] = strandReads[0][iBase];
                    n[1] = strandReads[1][iBase];
                    for (int jBase = 0; jBase < 4; jBase++) {
                        if (iBase != jBase) {
                            m[0] += strandReads[0][jBase];
                            m[1] += strandReads[1][jBase];
                        }
                    }
                    ratio[0] = 0.0;
                    ratio[1] = 0.0;
                    if (totStrand[0] > 0) {
                        ratio[0] = (1.0 * n[0]) / (1.0 * (n[0] + m[0]));
                    }
                    if (totStrand[1] > 0) {
                        ratio[1] = (1.0 * n[1]) / (1.0 * (n[1] + m[1]));
                    }
//                    System.out.println(n[0] + "\t" + m[0] + "\t" + n[1] + "\t" + m[1] + "\t" + Arrays.toString(ratio));
                    if (Math.max(ratio[0], ratio[1]) > topRatio) {
                        topRatio = Math.max(ratio[0], ratio[1]);
                        iConsensus = iBase;
                    }
                    double independent = - Math.log(n[0] + m[0] +1.0) - Math.log(n[1] + m[1] + 1.0);
                    double correlatedNum = CombinatoricsUtils.factorialLog(n[0]+m[0]) 
                            + CombinatoricsUtils.factorialLog(n[1]+m[1]) 
                            + CombinatoricsUtils.factorialLog(n[0]+n[1]) 
                            + CombinatoricsUtils.factorialLog(m[0]+m[1]);
                    double correlatedDen = CombinatoricsUtils.factorialLog(n[0]) 
                            + CombinatoricsUtils.factorialLog(n[1]) 
                            + CombinatoricsUtils.factorialLog(m[0]) 
                            + CombinatoricsUtils.factorialLog(m[1])
                            + CombinatoricsUtils.factorialLog(n[0] + n[1] + m[0] + m[1] + 1);
                    double correlated = correlatedNum - correlatedDen;
                    double posterior = 0.99 * Math.exp(correlated) / (0.99 * Math.exp(correlated) + 0.01 * Math.exp(independent));
                    minAmt[iBase] = posterior * inverseBeta(0.025, (0.5+reads[iBase]), (0.5+totReads-reads[iBase]));
                    if (posterior < 0.01) printThis = true;
                }
            }
            minAmt[iConsensus] = 1.0;
            if (printThis) {
                System.out.println(Arrays.toString(strandReads[0]) + "\t" + Arrays.toString(strandReads[1])
                + "\t" + Arrays.toString(minAmt));
            }
        }
    }
    
    boolean isConserved() {
        return conserved;
    }
    
    boolean hasData() {
        return hasData;
    }
    
    int getISite() {
        return iSite;
    }
    
    double[] getMinAmt() {
        return minAmt;
    }
    
    int[] getReads() {
        return reads;
    }
    
    public double inverseBeta(double p, double alpha, double beta){
        BetaDistribution bDist = new BetaDistribution(alpha, beta);
        double x = 0;
        double a = 0;
        double b = 1;
        double precision = 0.01; // converge until there is 2 decimal places precision

        while ((b - a) > precision) {
            x = (a + b) / 2;
            if (bDist.cumulativeProbability(x) > p) {
                b = x;
            } else {
                a = x;
            }
        }
        return x;
    }
    
}

