package rag.cluster;

import java.util.ArrayList;
import java.util.Arrays;


/**
 * Container for the reads at a single site
 *
 * @author rgoldst
 */
public class Site {
    private final GammaCalc gamma;
    int iSite;
    double[] estProbDiffBases = new double[5];
    int conservedBase = -9;
    boolean siteConserved = false;
    double totalLogLikelihood = 0.0;
    private int nTimePoints;
    private int nHaplo = 0;
    private int nBases = 0;
    private ArrayList<Assignment> assignmentVector;
    private ArrayList<Assignment> localAssignmentVector = new ArrayList<>();
    private double[] probAssignment = null;
    private double[] priorProb = new double[2];
    private int nAssignments = 0;
    private String[] baseString = {"A", "C", "G", "T"};
    private int[][][] strandReads = null; // [tp][strand][base] top two sets of reads on each strand
    private int[][] totStrand = null; // [tp][strand] number of reads on each strand
    private int[][] reads = null; // [tp][base]
    private int[] totReads = null; // [tp]
    private int[] timePointConservedBase = null;
    private boolean siteActive = false;
    private int nPresentBase = 0;   // number of present bases
    private boolean[] presentBase = new boolean[4];
    private boolean[] timePointConserved = null;
    private boolean[] timePointHasData = null;

    Site(int iSite, int nTimePoints, int nHaplo, ArrayList<Assignment> assignmentVector, GammaCalc gammaCalc) {
        this.gamma = gammaCalc;
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
                strandReads[iTimePoint][iStrand][iBase] = Integer.parseInt(words[(2 * iBase) + iStrand + 3]);
                reads[iTimePoint][iBase] += strandReads[iTimePoint][iStrand][iBase];
            }
//            strandReads[iTimePoint][0][iBase] = reads[iTimePoint][iBase]/2;
//            strandReads[iTimePoint][1][iBase] = reads[iTimePoint][iBase] - strandReads[iTimePoint][0][iBase];
            for (int iStrand = 0; iStrand < 2; iStrand++) {
                totStrand[iTimePoint][iStrand] += strandReads[iTimePoint][iStrand][iBase];
                totReads[iTimePoint] += strandReads[iTimePoint][iStrand][iBase];
            }


            if (reads[iTimePoint][iBase] > 0) {
                conservedBase = iBase;
                presentBase[iBase] = true;
            }
        }
        if (totReads[iTimePoint] > 0) {   // has data
            timePointHasData[iTimePoint] = true;
        }
        if (Math.max(Math.max(reads[iTimePoint][0], reads[iTimePoint][1]),
                Math.max(reads[iTimePoint][2], reads[iTimePoint][3])) == totReads[iTimePoint]) {  // timePointConserved site
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
        smellTest();
        return siteActive;
    }

    double assignHaplotypes(double[] alphaParams, double[] priors) {
        double alpha0 = alphaParams[0] * (1.0 - alphaParams[1]) / alphaParams[1];
        double alphaE = (1.0 - alphaParams[0]) * (1.0 - alphaParams[1]) / alphaParams[1];
        estProbDiffBases = new double[5];
        double logLikelihood = 0.0;
        if (siteConserved) {
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                for (int iStrand = 0; iStrand < 2; iStrand++) {
                    logLikelihood += priors[1] + this.gamma.logGamma(alpha0 + 3.0 * alphaE)
                            - this.gamma.logGamma(alpha0 + 3.0 * alphaE + totStrand[iTimePoint][iStrand])
                            + this.gamma.logGamma(alpha0 + totStrand[iTimePoint][iStrand])
                            - this.gamma.logGamma(alpha0);
                }
            }
            estProbDiffBases[1] = 1.0;
            return logLikelihood;
        }
        probAssignment = new double[localAssignmentVector.size()];
        double[] logLikelihoodAssign = new double[localAssignmentVector.size()];
        double sumProb = 0.0;
        int bestAssign = -999;
        double bestAssignVal = -1.0E20;
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            Assignment assignment = localAssignmentVector.get(iAssign);
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                logLikelihoodAssign[iAssign] += priors[localAssignmentVector.get(iAssign).nPresent]
                        + assignment.computeAssignmentLogLikelihood(iTimePoint, strandReads[iTimePoint],
                        reads[iTimePoint], totStrand[iTimePoint], siteConserved);
            }
            if (logLikelihoodAssign[iAssign] > bestAssignVal) {
                bestAssignVal = logLikelihoodAssign[iAssign];
                bestAssign = iAssign;
            }
        }


        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            probAssignment[iAssign] = Math.exp(logLikelihoodAssign[iAssign] - bestAssignVal);
            sumProb += probAssignment[iAssign];
            logLikelihood += probAssignment[iAssign];
            double nContrib = Math.exp(logLikelihoodAssign[iAssign] - bestAssignVal
                    - priors[localAssignmentVector.get(iAssign).nPresent]
                    + priors[localAssignmentVector.get(bestAssign).nPresent]);
            estProbDiffBases[localAssignmentVector.get(iAssign).nPresent]
                    += nContrib;
            estProbDiffBases[0] += nContrib;
        }
        logLikelihood = bestAssignVal + Math.log(logLikelihood);
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            probAssignment[iAssign] /= sumProb;
        }
        for (int nBase = 1; nBase < 5; nBase++) {
            estProbDiffBases[nBase] /= estProbDiffBases[0];
        }

        if (false) {
            System.out.print(iSite);
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                System.out.print("\t" + Arrays.toString(reads[iTimePoint]));
            }
            System.out.println("\t" + Arrays.toString(localAssignmentVector.get(bestAssign).assign) + "\t" + bestAssignVal);

        }
        return logLikelihood;
    }

    boolean smellTest() {
        double estProb = 0.0;
        double actProb = 0.0;
        boolean smellBad = false;
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            if (totReads[iTimePoint] * totStrand[iTimePoint][0] * totStrand[iTimePoint][1] > 0) {
                for (int iBase = 0; iBase < 4; iBase++) {
                    if (reads[iTimePoint][iBase] > 0) {
                        estProb += this.gamma.logGamma(reads[iTimePoint][iBase] + 0.5)
                                - this.gamma.logGamma(reads[iTimePoint][iBase] + 1.0) - 0.5723649;
                        actProb += -reads[iTimePoint][iBase] * 0.6931472 + this.gamma.logGamma(reads[iTimePoint][iBase] + 1)
                                - this.gamma.logGamma(1.0 + 0.5 * totReads[iTimePoint] * strandReads[iTimePoint][0][iBase] / totStrand[iTimePoint][0])
                                - this.gamma.logGamma(1.0 + 0.5 * totReads[iTimePoint] * strandReads[iTimePoint][1][iBase] / totStrand[iTimePoint][1]);
                    }
                }
            }
        }
        smellBad = (actProb - estProb < -11.512); // 11.512
        if (smellBad) {
            System.out.println("Rejected site\t" + iSite + "\t" + (actProb - estProb) + "\t" + smellBad);
        }
        return smellBad;
    }


    double computeSiteLogLikelihood(double[] alphaParams, double[] priors) {
        double alpha0 = alphaParams[0] * (1.0 - alphaParams[1]) / alphaParams[1];
        double alphaE = (1.0 - alphaParams[0]) * (1.0 - alphaParams[1]) / alphaParams[1];
        totalLogLikelihood = 0.0;
        if (siteConserved) {
            double g1 = priors[1] + this.gamma.logGamma(alpha0 + 3.0 * alphaE) - this.gamma.logGamma(alpha0);
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                int[] thisStrand = totStrand[iTimePoint];
                for (int iStrand = 0; iStrand < 2; iStrand++) {
                    totalLogLikelihood += g1
                            - this.gamma.logGamma(alpha0 + 3.0 * alphaE + thisStrand[iStrand])
                            + this.gamma.logGamma(alpha0 + thisStrand[iStrand]);
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
                if (probAssignment[iAssign] > 0.01) {
                    Assignment assignment = localAssignmentVector.get(iAssign);
                    logLikelihoodAssign[iAssign] += priors[localAssignmentVector.get(iAssign).nPresent]
                            + assignment.computeAssignmentLogLikelihood(iTimePoint, strandReads[iTimePoint],
                            reads[iTimePoint], totStrand[iTimePoint], siteConserved);
                    if (logLikelihoodAssign[iAssign] > bestAssignVal) {
                        bestAssignVal = logLikelihoodAssign[iAssign];
                        bestAssign = iAssign;
                    }
                }
            }
            for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
                if (probAssignment[iAssign] > 0.01) {
                    timePointLogLikelihood += probAssignment[iAssign] * Math.exp(logLikelihoodAssign[iAssign] - bestAssignVal);
                }
            }

            totalLogLikelihood += bestAssignVal + Math.log(timePointLogLikelihood);
        }
//        if (iSite == 230865) System.out.println("totalLogLikelihood B " + totalLogLikelihood);
        return totalLogLikelihood;
    }


    double computeSiteTimePointLogLikelihood(int iTimePoint, double[] alphaParams, double[] priors) {
        double alpha0 = alphaParams[0] * (1.0 - alphaParams[1]) / alphaParams[1];
        double alphaE = (1.0 - alphaParams[0]) * (1.0 - alphaParams[1]) / alphaParams[1];
        totalLogLikelihood = 0.0;
        if (siteConserved) {
            for (int iStrand = 0; iStrand < 2; iStrand++) {
                totalLogLikelihood += priorProb[0] + this.gamma.logGamma(alpha0 + 3.0 * alphaE)
                        - this.gamma.logGamma(alpha0 + 3.0 * alphaE + totStrand[iTimePoint][iStrand])
                        + this.gamma.logGamma(alpha0 + totStrand[iTimePoint][iStrand])
                        - this.gamma.logGamma(alpha0);
            }
            return totalLogLikelihood;
        }
        double[] logLikelihoodAssign = new double[localAssignmentVector.size()];
        int bestAssign = -999;
        double bestAssignVal = -1.0E20;
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            if (probAssignment[iAssign] > 0.01) {
                Assignment assignment = localAssignmentVector.get(iAssign);
                logLikelihoodAssign[iAssign] += priors[localAssignmentVector.get(iAssign).nPresent]
                        + assignment.computeAssignmentLogLikelihood(iTimePoint, strandReads[iTimePoint],
                        reads[iTimePoint], totStrand[iTimePoint], siteConserved);
                if (logLikelihoodAssign[iAssign] > bestAssignVal) {
                    bestAssignVal = logLikelihoodAssign[iAssign];
                    bestAssign = iAssign;
                }
            }
        }
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            if (probAssignment[iAssign] > 0.01) {
                totalLogLikelihood += probAssignment[iAssign] * Math.exp(logLikelihoodAssign[iAssign] - bestAssignVal);
            }
        }
        totalLogLikelihood = bestAssignVal + Math.log(totalLogLikelihood);
        return totalLogLikelihood;
    }


    double[][] getProbBase() {
        double[][] expectedFreq = new double[nHaplo][4];
        if (siteConserved) {
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                expectedFreq[iHaplo][conservedBase] = 1.0;
            }
        } else {
            for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
                Assignment assignment = localAssignmentVector.get(iAssign);
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    expectedFreq[iHaplo][assignment.assign[iHaplo]] += probAssignment[iAssign];
                }
            }
        }
        if (!siteConserved) {
            if (smellTest()) {
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    for (int iBase = 0; iBase < 4; iBase++) {
                        expectedFreq[iHaplo][iBase] = 0.0;
                    }
                }
            }
        }
        return expectedFreq;
    }


}
