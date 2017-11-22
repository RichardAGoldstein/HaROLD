/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Cluster;

import java.util.Arrays;
import org.apache.commons.math3.distribution.BetaDistribution;
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
    double[] minAmt = new double[4];   // minimum frequency of base
    
    boolean hasData = true;
    boolean conserved = false;
    
    Assembly(String inputLine){   // create data point and read in line of input file
        this.inputLine = inputLine;
        String[] words = inputLine.split(",");
        iSite = Integer.parseInt(words[0]);

        for (int iBase = 0; iBase < 4; iBase++) {   // compute various sums of reads
            for (int iStrand = 0; iStrand < 2; iStrand++) {
                strandReads[iStrand][iBase]=Integer.parseInt(words[(2*iBase)+iStrand+3]);
                totStrand[iStrand] += strandReads[iStrand][iBase];
                reads[iBase] += strandReads[iStrand][iBase];
                totReads += strandReads[iStrand][iBase];
            }
        }
        if (totReads == 0) {   // no data
            hasData = false;
        } else if (Math.max(Math.max(reads[0],reads[1]),Math.max(reads[2],reads[3])) == totReads) {  // conserved site
            conserved = true;
            int iConsensus = -999;
            for (int iBase = 0; iBase < 4; iBase++) {
                if (reads[iBase] == totReads) {
                    iConsensus = iBase;
                }
            }
            if (iConsensus < 0) {
                System.out.println("Something wrong in Assembly 1");
                System.exit(1);
            }
            minAmt[iConsensus] = 1.0 + 1.0E-4 * Cluster.random.nextDouble();
        } else {    // compute relative amounts of 4 different bases
            double topRatio = 0.0;
            int iConsensus = -999;
            double[] posterior = new double[4];
            for (int iBase = 0; iBase < 4; iBase++) {
                if (reads[iBase] > 0) {
                    int[] n = new int[2];
                    int[] m = new int[2];
                    double[] ratio = new double[2];
                    n[0] = strandReads[0][iBase];
                    n[1] = strandReads[1][iBase];
                    m[0] = totStrand[0] - n[0];
                    m[1] = totStrand[1] - n[1];
                    ratio[0] = 0.0;
                    ratio[1] = 0.0;
                    if (totStrand[0] > 0) {
                        ratio[0] = (1.0 * n[0]) / (1.0 * (n[0] + m[0]));   // ratio on two strands
                    }
                    if (totStrand[1] > 0) {
                        ratio[1] = (1.0 * n[1]) / (1.0 * (n[1] + m[1]));
                    }
                    if (Math.max(ratio[0], ratio[1]) > topRatio) {   // Assuming base with maximum frequency on strand is consensus 
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
                    posterior[iBase] = 0.99 * Math.exp(correlated) / (0.99 * Math.exp(correlated) + 0.01 * Math.exp(independent));
                }
            }
            posterior[iConsensus] = 1.0;
            for (int iBase = 0; iBase < 4; iBase++) {
                if (reads[iBase] > 0) {
                    minAmt[iBase] = posterior[iBase] * inverseBeta(0.025, (0.5+reads[iBase]), (0.5+totReads-reads[iBase]));
                    minAmt[iBase] *= (1.0 + 1.0E-4 * Cluster.random.nextDouble());   // Add a small amount so no numbers are exactly equal
                }
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
    
    int[][] getStrandReads() {
        return strandReads;
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

