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
public class DataSet implements MultivariateFunction{
    
    Hashtable<Integer, Site> siteHash = new Hashtable<Integer, Site>();  // Data of sites labeled by site number
    Vector<Integer> activeSiteVector = new Vector<Integer>();  // List of sites that are actively considered
    Vector<Integer> allSiteVector = new Vector<Integer>();// List of all sites  
    int nHaplo = 3; // Number of haplotypes
    Vector<Assignment> assignmentVector = null;
    boolean addFlat = false;  // Add a 'garbage' model for random outliers  ***NOT IMPLEMENTED***
    int nTimePoints = 0;   // Number of time points
    
    int iter = 0;  // Initialise count of iterations   
    int iStage = 0;
    int iTimePoint = 0;
    
    int optType = 0;
    int optTimePoint = 0;
    
    double[][] currentAlphaHap = null;
    double currentAlpha_C = 1.0;
    double currentAlpha_E = 0.01;
    
    int iCount = 0;

    DataSet(String fileNameFile, int nHaplo, Vector<Assignment> assignmentVector, boolean addFlat) {  // Read in data
        this.nHaplo = nHaplo;
        this.assignmentVector = assignmentVector;
        this.addFlat = addFlat;
        if (addFlat)  {
            System.out.println("addFlat not yet implemented");
            System.exit(1);
        }
        Vector<String> fileNameVector = new Vector<String>(); // list of files to be read in one for each time point
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
                activeSiteVector.add(iSite);
            } 
        }
        currentAlphaHap = new double[nTimePoints][nHaplo];
        for (int iTimePoint = 0; iTimePoint< nTimePoints; iTimePoint++) {
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                currentAlphaHap[iTimePoint][iHaplo] = iHaplo + 1.0;
            }
        }
        currentAlpha_C = 1.0;
        currentAlpha_E = 0.001;
    }
    
    double computeTotalLogLikelihood(double[][] alphaHap, double alpha_C, double alpha_E) {
        for (Assignment assignment : assignmentVector) {
            assignment.setAlphas(alphaHap, alpha_C, alpha_E);
        }
        double totalLogLikelihood = 0.0;
        for (int iSite : activeSiteVector) {
            Site site = siteHash.get(iSite);
            totalLogLikelihood += site.computeSiteLogLikelihood();
        }
        return totalLogLikelihood;
    }
    
    void setOptType(int optType, int optTimePoint) {
        this.optType = optType;
        this.optTimePoint = optTimePoint;
        System.out.println(optType + "\t" + optTimePoint);
        iCount = 0;
    }
    
    
    
    
    int getNTimePoints() {
        return nTimePoints;
    }
//    
//    void setStage(int iStage, double alpha_e, double beta_e, double S, double F0, double[][] alpha) {
//        iter = 0;
//        this.iStage = iStage;
//        for (int iSite : activeSiteVector) {
//            siteHash.get(iSite).setParams(iStage, alpha_e, beta_e, S, F0, alpha);   
//        }
//    }
//    
//    void setITimePoint(int iTimePoint) {
//        this.iTimePoint = iTimePoint;
//    }
//    
    public double value(double[] params) {
        
        if (optType == 0) {
            double val = computeTotalLogLikelihood(currentAlphaHap, params[0], params[1]);
            if (iCount % 100 == 0) {
                System.out.print("xxx\t" + Arrays.toString(params));
                System.out.println("\t" + val);
            }
            iCount++;
            return -val;
        } else if (optType == 1) {
            double[][] newAlphaHap = new double[nTimePoints][nHaplo];
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    if (iTimePoint != optTimePoint){
                        newAlphaHap[iTimePoint][iHaplo] = currentAlphaHap[iTimePoint][iHaplo];
                    } else {
                        newAlphaHap[iTimePoint][iHaplo] = params[iHaplo];
                    }
                }
            }
            double val = computeTotalLogLikelihood(newAlphaHap, currentAlpha_C, currentAlpha_E);
            if (iCount % 100 == 0) {
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
}    
    


