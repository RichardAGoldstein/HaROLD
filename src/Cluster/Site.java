/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Cluster;

import java.util.Arrays;
import java.util.Vector;
import org.apache.commons.math3.special.Gamma;


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
    Vector<Assignment> assignmentVector = null;
    Vector<Assignment> localAssignmentVector = new Vector<>();
    double[] probAssignment = null;
    int nAssignments = 0;
    String[] baseString = {"A", "C", "G", "T"};
    
    
    int[][][] strandReads = null; // [tp][strand][base] top two sets of reads on each strand
    int[][] totStrand = null; // [tp][strand] number of reads on each strand
    int[][] reads = null; // [tp][base]
    int[] totReads = null; // [tp]
    int[] timePointConservedBase = null; 

    
    boolean siteActive = false;
    int nPresentBase = 0;   // number of present bases
    boolean[] presentBase = new boolean[4];
    boolean siteConserved = false;
    boolean[] timePointConserved = null;
    boolean[] timePointHasData = null;
    
    Site(int iSite, int nTimePoints, int nHaplo, Vector<Assignment> assignmentVector) {
        this.iSite = iSite;
        this.nTimePoints = nTimePoints;
        this.nHaplo = nHaplo;
        this.assignmentVector = assignmentVector;
        nAssignments = assignmentVector.size();
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

    
    boolean isActive() {
        siteActive = false;
        for (int iBase = 0; iBase < 4; iBase++) {
            if (presentBase[iBase]) {
                nPresentBase++;
                siteActive = true;
            }
        }
        siteConserved = (nPresentBase == 1);
        for (Assignment assignment : assignmentVector) {
            boolean addThis = true;
            for (int iBase = 0; iBase < 4; iBase++) {
                if (assignment.presentBase[iBase] && !presentBase[iBase]) {
                    addThis = false;
                }
            }
            if (addThis) {
                localAssignmentVector.add(assignment);
            }
        }
        return siteActive;
    }
    
    double assignHaplotypes() {
        double logLikelihood = 0.0;
        probAssignment = new double[localAssignmentVector.size()];
        double[] logLikelihoodAssign = new double[localAssignmentVector.size()];
        double sumProb = 0.0;
        int bestAssign = -999;
        double bestAssignVal = -1.0E20;
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            Assignment assignment = localAssignmentVector.get(iAssign);
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                logLikelihoodAssign[iAssign] += 
                        assignment.computeAssignmentLogLikelihood(iTimePoint, strandReads[iTimePoint], 
                                reads[iTimePoint], totStrand[iTimePoint], siteConserved);
            }
            if (logLikelihoodAssign[iAssign] > bestAssignVal) {
                bestAssignVal = logLikelihoodAssign[iAssign];
                bestAssign = iAssign;
            }
        }
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            probAssignment[iAssign] = Math.exp(logLikelihoodAssign[iAssign]-bestAssignVal);
            sumProb += probAssignment[iAssign];
            logLikelihood += probAssignment[iAssign];
        }
        logLikelihood = bestAssignVal + Math.log(logLikelihood);
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            probAssignment[iAssign] /= sumProb;
        }
        if (false && Cluster.verbose) {
            System.out.print(iSite);
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                System.out.print("\t" + Arrays.toString(reads[iTimePoint]));
            }
            System.out.println("\t" + Arrays.toString(localAssignmentVector.get(bestAssign).assign) + "\t" + bestAssignVal);
            
        }
        return logLikelihood;
    }
    
    
    double computeSiteLogLikelihood(double alpha0, double alphaE) {
        double totalLogLikelihood = 0.0; 
        if (siteConserved) {
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                for (int iStrand = 0; iStrand < 2; iStrand++) {
                    totalLogLikelihood = Gamma.logGamma(alpha0 + 3.0 * alphaE)
                            - Gamma.logGamma(alpha0 + 3.0 * alphaE + totStrand[iTimePoint][iStrand])
                            + Gamma.logGamma(alpha0 + totStrand[iTimePoint][iStrand])
                            - Gamma.logGamma(alpha0);
                }
            }
            return totalLogLikelihood;
        }

        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            double[] logLikelihoodAssign = new double[localAssignmentVector.size()];
            double timePointLogLikelihood = 0.0;
            int bestAssign = -999;
            double bestAssignVal = -1.0E20;
            for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
                if (probAssignment[iAssign] > 0.001) {
                    Assignment assignment = localAssignmentVector.get(iAssign);
                    logLikelihoodAssign[iAssign] += 
                            assignment.computeAssignmentLogLikelihood(iTimePoint, strandReads[iTimePoint], 
                                    reads[iTimePoint], totStrand[iTimePoint], siteConserved);
                    if (logLikelihoodAssign[iAssign] > bestAssignVal) {
                        bestAssignVal = logLikelihoodAssign[iAssign];
                        bestAssign = iAssign;
                    }
                }
            }
            for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
                if (probAssignment[iAssign] > 0.001) {
                    timePointLogLikelihood += probAssignment[iAssign] * Math.exp(logLikelihoodAssign[iAssign]-bestAssignVal);
                }
            }
            totalLogLikelihood += bestAssignVal + Math.log(timePointLogLikelihood);
        }
        return totalLogLikelihood;
    }
    
    
    double computeSiteTimePointLogLikelihood(int iTimePoint, double alpha0, double alphaE) { 
        double totalLogLikelihood = 0.0;
        if (siteConserved) {
            for (int iStrand = 0; iStrand < 2; iStrand++) {
                totalLogLikelihood = Gamma.logGamma(alpha0 + 3.0 * alphaE)
                        - Gamma.logGamma(alpha0 + 3.0 * alphaE + totStrand[iTimePoint][iStrand])
                        + Gamma.logGamma(alpha0 + totStrand[iTimePoint][iStrand])
                        - Gamma.logGamma(alpha0);
            }
            return totalLogLikelihood;
        }
        double[] logLikelihoodAssign = new double[localAssignmentVector.size()];
        int bestAssign = -999;
        double bestAssignVal = -1.0E20;
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            if (probAssignment[iAssign] > 0.001) {
                Assignment assignment = localAssignmentVector.get(iAssign);
                    logLikelihoodAssign[iAssign] += 
                            assignment.computeAssignmentLogLikelihood(iTimePoint, strandReads[iTimePoint], 
                                    reads[iTimePoint], totStrand[iTimePoint], siteConserved);
                if (logLikelihoodAssign[iAssign] > bestAssignVal) {
                    bestAssignVal = logLikelihoodAssign[iAssign];
                    bestAssign = iAssign;
                }
            }
        }
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            if (probAssignment[iAssign] > 0.001) {
                totalLogLikelihood += probAssignment[iAssign] * Math.exp(logLikelihoodAssign[iAssign]-bestAssignVal);
            }
        }
        totalLogLikelihood = bestAssignVal + Math.log(totalLogLikelihood);
        return totalLogLikelihood;
    }

}
