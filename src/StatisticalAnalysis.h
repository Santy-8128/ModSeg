
#ifndef __STATISTICAL_ANALYSIS_H
#define __STATISTICAL_ANALYSIS_H

#include <ctime>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "StringBasics.h"
#include "Pedigree.h"
#include "PedigreeData.h"
using namespace std;



class StatAnalysis
{
    protected:


    public:
        String mapFile;
        String pedFile;
        String condPedFile;
        int condition,model,parentEffect,phenoStrat,sexStrat;
        vector<double> MLE;
        vector<string> ParamName;
        double MaxLogLikelihood;
        double base_constant;
        double pMin,pMax,pMid;
        double critic;
        vector<double> baselineHazard;
        vector<double> cummulativeBaselineHazard;
        vector<double> cummulativeRisk;

        vector<double> mlePvalue;

        void printCumRisk();
        void printCumHazard();

        void InvInformation();
        vector<vector<double> > InformationMatrix;
        vector<vector<double> > VarianceMLE;
        bool singular;

        StatAnalysis(String pedfile,String mapfile,String condPedfile,map<string,double> &Param);


        void PerformAnalysis();
        void resultPrint();
        //friend class Likelihood;
        friend class Node;
};






#endif
