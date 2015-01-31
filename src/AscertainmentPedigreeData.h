


#ifndef __ASCERTAINMENT_PEDIGREE_DATA_H
#define __ASCERTAINMENT_PEDIGREE_DATA_H

#include <ctime>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "Pedigree.h"
#include <string>
#include <vector>
#include <map>
using namespace std;



class AscertainmentPedigreeData
{
    public:

        int DataModel;

        // 0 = for case-control study, No Age Strat, No Sex Strat
        // 1 = for survival analysis study, No Age Strat, No Sex Strat
        // 0 = for case-control study, No Age Strat, No Sex Strat
        // 0 = for case-control study, No Age Strat, No Sex Strat
        // 0 = for case-control study, No Age Strat, No Sex Strat
        // 0 = for case-control study, No Age Strat, No Sex Strat
        // 0 = for case-control study, No Age Strat, No Sex Strat
        // 0 = for case-control study, No Age Strat, No Sex Strat

        string MarkerName;
        string Allele1,Allele2;
        string PhenoName;
        double alleleFreq1,baseline;
        int censorIndicator;
        double phenoMin,phenoMax,phenoMid;;
        vector<vector<string> > peddata;
        vector<vector<string> > Condpeddata;


        int NoOfPheno;
        int TotalSample1,TotalSample2;
        int NoOfPedigree1,NoOfPedigree2;
        vector<string> PedigreeId1,PedigreeId2;



        vector<Pedigree> Data;
        vector<Pedigree> CondData;
        AscertainmentPedigreeData(char* pedfile,char* condpedfile,char* mapfile,int model,double baselineHazard);
        AscertainmentPedigreeData()
        {

        };
        void ImportData();

        void readPedigreeFile(char* file,int cond);
        double operator[](double x);
        void printMapFileSummary();
        void readMapFile(char* mapfile);
        double operator()(vector<double> x);
        //double FindMLE();


    //friend class Likelihood;
    friend class Node;
};






#endif
