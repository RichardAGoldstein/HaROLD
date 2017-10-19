/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Cluster;

import java.util.*;
import org.apache.commons.math3.distribution.BetaDistribution;
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
            }
            reads[iBase] = strandReads[0][iBase]+strandReads[1][iBase];
            totReads += reads[iBase];
        }
        if (totReads == 0) {
            hasData = false;
        } else if (Math.max(Math.max(reads[0],reads[1]),Math.max(reads[2],reads[3])) == totReads) {
            conserved = true;
        } // estimate lower bound on number of reads of each base
        for (int iBase = 0; iBase < 4; iBase++) {
            if (reads[iBase] > 0) {   // downweight for excessive strand bias
                double ratio = (strandReads[0][iBase]+1.0)/(reads[iBase]+2.0);
                double zScore = Math.abs(ratio-0.5)*Math.sqrt(2.0*reads[iBase]);
                double pVal = (1.0 - Erf.erf(zScore));
//                System.out.println(iSite + "\t" + Arrays.toString(reads) + "\t" + pVal);
                minAmt[iBase] = (pVal/(0.01 + pVal))  // consider lowest CI value
                    * inverseBeta(0.025, (0.5+reads[iBase]), (0.5+totReads-reads[iBase]));
            }
        }
        System.out.println(iSite + "\t" + Arrays.toString(reads) + "\t" + Arrays.toString(minAmt));
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

