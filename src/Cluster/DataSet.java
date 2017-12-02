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
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;

/**
 *
 * @author rgoldst
 */
public class DataSet implements MultivariateFunction, UnivariateFunction {
    Vector<Site> activeSiteVector = new Vector<>();  // List of sites that are actively considered
    Vector<Site> variableSiteVector = new Vector<>(); // List of all variable sites
    int nHaplo = 3; // Number of haplotypes
    Vector<Assignment> assignmentVector = null;
    int nTimePoints = 0;   // Number of time points
    double[] alphaParams = new double[2];
         
    int optType = 0;
    int optTimePoint = 0;    
    int iCount = 0;

    DataSet(String fileNameFile, int nHaplo, Vector<Assignment> assignmentVector) {  // Read in data
        this.nHaplo = nHaplo;
        this.assignmentVector = assignmentVector;
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
                if (!site.siteConserved) {
                    variableSiteVector.add(site);
                }
            } 
        }
    }
    
    double computeTotalLogLikelihood(double[] alphaParams) {
        double alpha0 = alphaParams[0];
        double alphaE = alphaParams[1];
        double totalLogLikelihood = 0.0;
        if (optType == 0) {
            for (Site site : activeSiteVector) {
                totalLogLikelihood += site.computeSiteLogLikelihood(alpha0, alphaE);
            }
            return totalLogLikelihood;
        } else {
            for (Site site : variableSiteVector) {
                totalLogLikelihood += site.computeSiteTimePointLogLikelihood(optTimePoint, alpha0, alphaE);
            }
            return totalLogLikelihood;            
        }
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
 
    void assignHaplotypes(double[][] piParams, double[] alphaParams) {
        setParams(piParams, alphaParams);
        double totalLogLikelihood = 0.0;
        for (Site site : variableSiteVector) {
            totalLogLikelihood += site.assignHaplotypes();
        }
        System.out.println("zzz\t" + totalLogLikelihood);
    }
    
    void setParams(double[][] piParams, double[] alphaParams) {
        double[][] piHap = computePiHap(piParams);
        for (Assignment assignment : assignmentVector) {
            assignment.setParams(piHap, alphaParams);
        }
    }
    
    double[][] computePiHap(double[][] piParams) {
        double[][] piHap = new double[nTimePoints][4];
        double[] remaining = new double[nTimePoints];
        Arrays.fill(remaining, 1.0);
        for (int iHaplo = 0; iHaplo < nHaplo-1; iHaplo++) {
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                piHap[iTimePoint][iHaplo] = remaining[iTimePoint] * piParams[iTimePoint][iHaplo];
                remaining[iTimePoint] -= piHap[iTimePoint][iHaplo];
            }
        }
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            piHap[iTimePoint][nHaplo-1] = remaining[iTimePoint];
        }
        return piHap;
    }
    
    double[] computePiHap(double[] piParams) {
        double[] piHap = new double[4];
        double remaining = 1.0;
        for (int iHaplo = 0; iHaplo < nHaplo-1; iHaplo++) {
            piHap[iHaplo] = remaining * piParams[iHaplo];
            remaining -= piHap[iHaplo];
        }
        piHap[nHaplo-1] = remaining;
        return piHap;
    }
    
    
    void setAlphaParams(double[] alphaParams) {
        this.alphaParams = alphaParams;
        for (Assignment assignment : assignmentVector) {
            assignment.setAlphas(alphaParams);
        }
    }
    
    void setPiHap(int iTimePoint, double[] piParams) {
        double[] piHap = computePiHap(piParams);
        for (Assignment assignment : assignmentVector) {
            assignment.setPiHap(iTimePoint, piHap);
        }
    }
    
    public double value(double[] params) {     
        if (optType == 0) {
            setAlphaParams(params);
            double val = computeTotalLogLikelihood(params);
            if (iCount % 10 == 0) {
                System.out.print("yyy\t" + Arrays.toString(params));
                System.out.println("\t" + val);
            }
            iCount++;
            return -val;
        } else if (optType == 1) {
            setPiHap(optTimePoint, params);
            double val = computeTotalLogLikelihood(alphaParams);
            if (iCount % 10 == 0) {
                System.out.print("xxx\t" + Arrays.toString(alphaParams) + "\t" + Arrays.toString(params));
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
    


