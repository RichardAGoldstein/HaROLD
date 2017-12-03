/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Cluster;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Vector;
import org.apache.commons.math3.analysis.MultivariateFunction;


/**
 *
 * @author rgoldst
 */
public class DataSet implements MultivariateFunction {
    Vector<Site> activeSiteVector = new Vector<>();  // List of sites that are actively considered
    Vector<Site> variableSiteVector = new Vector<>(); // List of all variable sites
    Vector<Site> reducedSiteVector = new Vector<>();  
    int nHaplo = 3; // Number of haplotypes
    Vector<Assignment> assignmentVector = null;   // Vectir if assignments
    int nTimePoints = 0;   // Number of time points
    double[] currentAlphaParams = new double[2];   // alpha0 and alphaE
    double[][] currentPiHap = null;
    double useFrac = 1.0;
    boolean firstRound = true;
         
    String[] baseString = {"A", "C", "G", "T"};
        
    int optType = 0;  // 0 for optimising alpha0 and alphaE, 1 for optimising haplotype frequencies
    int optTimePoint = 0;   // if optType = 1, what timePoint is being optimised 
    int iCount = 0;  // How many iterations of optimiser have been finished
    
    double currentLogLikelihood = 0.0;

    DataSet(String fileNameFile, int nHaplo, Vector<Assignment> assignmentVector, double useFrac) {  // Read in data
        this.nHaplo = nHaplo;
        this.assignmentVector = assignmentVector;
        this.useFrac = useFrac;
        Vector<String> fileNameVector = new Vector<String>(); // list of files to be read in one for each time point
        Vector<Integer> allSiteVector = new Vector<>();// List of all sites 
        Hashtable<Integer, Site> siteHash = new Hashtable<Integer, Site>();  // Data of sites labeled by site number
        
        try {
            FileReader file = new FileReader(fileNameFile); 
            BufferedReader buff = new BufferedReader(file);;
            boolean eof = false;
            while (!eof) {
                String line = buff.readLine();
                if (line == null) {
                    eof = true;
                } else {
                    fileNameVector.add(line); 
                }
            }
        }
        catch (IOException e) {
            System.out.println("Error: File not found (IO error)");
            System.exit(1);
        }
        nTimePoints = fileNameVector.size();  // Number of timepoints = number of files

        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {  // read in data files
            String dataFile = fileNameVector.get(iTimePoint);
            try{
                FileReader file = new FileReader(dataFile);
                BufferedReader buff = new BufferedReader(file);
                String line = "";
                boolean eof = false;
                while (!eof) {
                    line = buff.readLine();
                    if (line == null) {
                            eof = true;
                    } else {
                        if (line.contains("Position")) {
                            line = buff.readLine();
                        }
                        int iSite = Integer.parseInt(line.split(",")[0]);
                        if (!allSiteVector.contains(iSite)) {   // list of sites that contain data
                            allSiteVector.add(iSite);
                            Site newSite = new Site(iSite, nTimePoints, nHaplo, assignmentVector); // create new site if needed
                            siteHash.put(iSite, newSite);
                        }
                        siteHash.get(iSite).addTimePoint(iTimePoint, line);  // add datapoint to site
                    }
                }
            }     
            catch (IOException e) {
                System.out.println("Error: File not found (IO error)");
                System.exit(1);
            }
        }
        
        for (int iSite : allSiteVector) {  // Create activeSiteVector
            Site site = siteHash.get(iSite);
            if (site.isActive()) {    // do simple sums
                activeSiteVector.add(site);
                if (Cluster.random.nextDouble() < useFrac) {
                    reducedSiteVector.add(site);
                }
                if (!site.siteConserved) {
                    variableSiteVector.add(site);
                }
            } 
        }
    }
    
    double computeTotalLogLikelihood() {
        double totalLogLikelihood = 0.0;
        if (optType == 0 && (!firstRound || useFrac > 0.99999)) {
            for (Site site : activeSiteVector) {
                totalLogLikelihood += site.assignHaplotypes(currentAlphaParams);
            }
            return totalLogLikelihood;
        } else if (optType == 0) {
            for (Site site : reducedSiteVector) {
                totalLogLikelihood += site.assignHaplotypes(currentAlphaParams);
            }
            return totalLogLikelihood;
        } else {
            for (Site site : variableSiteVector) {
                totalLogLikelihood += site.computeSiteTimePointLogLikelihood(optTimePoint, currentAlphaParams);
            }
            return totalLogLikelihood;            
        }
    } 
    
    void setOptType(int optType, int optTimePoint, double[][] hapParams, double[] alphaParams, boolean firstRound) {
        this.optType = optType;
        this.optTimePoint = optTimePoint;
        this.firstRound = firstRound;
        updateAllParams(hapParams, alphaParams);
        System.out.println(optType + "\t" + optTimePoint);
        iCount = 0;
    }
 
    double assignHaplotypes() {
        if (Cluster.verbose) {
            System.out.print(Arrays.toString(currentAlphaParams));
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                System.out.print("  " + Arrays.toString(currentPiHap[iTimePoint]));
            }
            System.out.println();
        }
        currentLogLikelihood = 0.0;
        for (Site site : activeSiteVector) {
            currentLogLikelihood += site.assignHaplotypes(currentAlphaParams);
        }
        System.out.println("zzz\t" + currentLogLikelihood);
        return currentLogLikelihood;
    }

    /**
     *
     * Update to new values of hapParams and alphaParams
     */    
    void updateAllParams(double[][] hapParams, double[] alphaParams) {
        currentAlphaParams = alphaParams;
        currentPiHap = computePiHap(hapParams);
        for (Assignment assignment : assignmentVector) {
            assignment.setAllParams(currentPiHap, currentAlphaParams);
        }
    }
    
    /**
     *
     * Update to new values of alphaParams
     */  
    void updateAlphaParams(double[] alphaParams) {
        currentAlphaParams = alphaParams;
        for (Assignment assignment : assignmentVector) {
            assignment.setAllParams(currentPiHap, currentAlphaParams);
        }
    }

    /**
     *
     * Update to new values of hapParams for single timepoint
     */  
    void updateSingleHapParams(int iTimePoint, double[] hapParams) {
        currentPiHap[iTimePoint] = computePiHap(hapParams);
        for (Assignment assignment : assignmentVector) {
            assignment.setSinglePiHap(iTimePoint, currentPiHap[iTimePoint]);
        }
    }
    
    public double value(double[] params) {     
        if (optType == 0) {
            updateAlphaParams(params);
            double val = computeTotalLogLikelihood();
            if (iCount % 10 == 0) {
                if (Cluster.verbose) {
                    System.out.print(Arrays.toString(currentAlphaParams));
                    for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                            System.out.print("  " + Arrays.toString(currentPiHap[iTimePoint]));
                    }
                    System.out.println();
                }
                System.out.print("yyy\t" + Arrays.toString(params));
                System.out.println("\t" + val);
            }
            iCount++;
            return -val;
        } else if (optType == 1) {
            updateSingleHapParams(optTimePoint, params);
            double val = computeTotalLogLikelihood();
            if (iCount % 10 == 0) {
                if (Cluster.verbose) {
                    System.out.print(Arrays.toString(currentAlphaParams));
                    for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                            System.out.print("  " + Arrays.toString(currentPiHap[iTimePoint]));
                    }
                    System.out.println();
                }
                System.out.print("xxx\t" + Arrays.toString(params));
                System.out.println("\t" + val);
            }
            iCount++;
            return -val;
        }
        System.out.println("Error in optimisation");
        System.exit(1);
        return 0.0;
    }
    
    public double value(double singleParam) {  
        double[] params = new double[1];
        params[0] = singleParam;
        return value(params);
    }
  
    void printResults() {
        System.out.println("\n\n\nResults for nHaplotypes = " + nHaplo);
        int nParams = 2 + (nHaplo-1)*nTimePoints;
        System.out.println("Number of adjustable parameters: " + nParams);
        System.out.println("Final likelihood: " + currentLogLikelihood);
        System.out.println("Dirichlet parameters for errors: " + currentAlphaParams[0] + "\t" + currentAlphaParams[1]);
        System.out.println("Error rate: " + (currentAlphaParams[1]/(currentAlphaParams[0]+currentAlphaParams[1])));
        System.out.println("Haplotype frequencies");
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            System.out.print(iTimePoint);
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                System.out.print("\t" + currentPiHap[iTimePoint][iHaplo]);          
            }
            System.out.println();
        }
        System.out.println("Haplotypes");
        int nSites = activeSiteVector.size();
        int[][] bestBase = new int[nSites][nHaplo];
        double[][] probBestBase = new double[nSites][nHaplo];
        for (int cSite = 0; cSite < nSites; cSite++) {
            Site site = activeSiteVector.get(cSite);
            double[][] probBase = site.getProbBase();
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                for (int iBase = 0; iBase < 4; iBase++) {
                    if (probBase[iHaplo][iBase] > probBestBase[cSite][iHaplo]) {
                        probBestBase[cSite][iHaplo] = probBase[iHaplo][iBase];
                        bestBase[cSite][iHaplo] = iBase;
                    }
                }
            }
        }
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
            System.out.print("Haplo_" + iHaplo + "\t");
            for (int cSite = 0; cSite < nSites; cSite++) {
                System.out.print(baseString[bestBase[cSite][iHaplo]]);              
            }
            System.out.println();
        }
        for (int cSite = 0; cSite < nSites; cSite++) {
            Site site = activeSiteVector.get(cSite);
            System.out.print(site.iSite);
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                System.out.print("\t" + Arrays.toString(site.reads[iTimePoint]));
            }
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                System.out.print("\t" + baseString[bestBase[cSite][iHaplo]]);              
            }
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                System.out.format("\t%6.4f", probBestBase[cSite][iHaplo]);              
            }
            System.out.println();
        }
    }
    
//    void printHaplotypes(double alpha_e, double beta_e, double S, double F0) {
//        System.out.println("zzz\nzzz\nzzz");
//        for (int iSite : activeSiteVector) {
//            Site site = siteHash.get(iSite);
//            site.computeLogLikelihood(alpha_e, beta_e, S, F0, true);   
//        }
//    }
//    
//    
//    void reAdjust(double alpha_e, double beta_e, double S, double F0) {
//        for (int iSite : activeSiteVector) {
//            Site site = siteHash.get(iSite);
//            site.computeLogLikelihood(alpha_e, beta_e, S, F0, false);   
//        }
//    }
//    
//    public double alphaValue(double[] params) {
//        double value = 0.0;
//        for (int iSite : activeSiteVector) {
//            Site site = siteHash.get(iSite);
//            value += site.computeLogLikelihoodAlpha(params, iTimePoint);   
//        }
//        return -value;
//    }
//
//    
//    public double globalValue(double[] params) {
//        double value = 0.0;
//        // Distribute parameters as appropriately
//        double alpha_e = params[0];
//        double S = params[1];
//        double F0 = params[2];
//        double beta_e = 10.0;
//        for (int iSite : activeSiteVector) {
//            Site site = siteHash.get(iSite);
//            value += site.computeLogLikelihood(alpha_e, beta_e, S, F0, false);   
//        }
//        if (iter%10 == 0) {  // Print stuff
//            System.out.format("zzz\t%d\t%10.4f\t%12.6f  %12.6f\t%10.8f\t%10.8f      ",
//                    iter, value, alpha_e, beta_e, S, F0);
//            System.out.println();
//        }
//        iter++;
//        return -value;
//    }
    
        /**
     *
     * Compute new values of piHap for all time points
     */    
    double[][] computePiHap(double[][] hapParams) {
        double[][] piHap = new double[nTimePoints][nHaplo];
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            piHap[iTimePoint] = computePiHap(hapParams[iTimePoint]);
        }
        return piHap;
    }
    
    /**
     *
     * compute new values of piHap for single time point
     */    
    double[] computePiHap(double[] hapParams) {
        double[] piHap = new double[nHaplo];
        double remaining = 1.0;
        for (int iHaplo = 0; iHaplo < nHaplo-1; iHaplo++) {
            piHap[iHaplo] = remaining * hapParams[iHaplo];
            remaining -= piHap[iHaplo];
        }
        piHap[nHaplo-1] = remaining;
        return piHap;
    }

    int getNTimePoints() {
        return nTimePoints;
    }
    
}    
    


