/*
 *  Copyright (C) 2013  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <sstream>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include "Likelihood.h"
#include "Pedigree.h"
#include "PedigreeData.h"
#include "Parameters.h"
#include "Optimization.h"
#include "Numerical.h"
#include "StringBasics.h"
# include "StatisticalAnalysis.h"

using namespace boost::math;
using namespace std;



void vcfVersion();
void description();
void usage();

int main(int argc, char ** argv)
{
    String pedfile = "";
    String mapfile = "";
    String condpedfile = "";
    int model=1,parentEffect=0,condition=0,sexStrat=0,phenoStrat=0;
    String outfile = "ModSeg.out";
    double baseHazard=0.0001;




    bool params = false;



    // Read in the parameters.
    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("ped", &pedfile)
        LONG_STRINGPARAMETER("map", &mapfile)
        LONG_INTPARAMETER("model", &model)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_DOUBLEPARAMETER("baseHazard", &baseHazard)
        LONG_INTPARAMETER("phenoStrat", &phenoStrat)
        LONG_INTPARAMETER("condition", &condition)
        LONG_INTPARAMETER("sexStrat", &sexStrat)
        LONG_INTPARAMETER("parentEffect", &parentEffect)
        LONG_STRINGPARAMETER("condped", &condpedfile)
        LONG_STRINGPARAMETER("out", &outfile)
        LONG_PARAMETER("params", &params)
        END_LONG_PARAMETERS();

    inputParameters.Add(new LongParameters ("Input Parameters",
                                            longParameterList));

    inputParameters.Read(argc, &(argv[0]));




    // Check that all files were specified.
    if(pedfile == "")
    {
        usage();
        inputParameters.Status();
        std::cerr << "Missing \"--ped\", a required parameter.\n\n";
        return(-1);
    }


    if(mapfile == "")
    {
        usage();
        inputParameters.Status();
        std::cerr << "Missing \"--map\", a required parameter.\n\n";
        return(-1);
    }


    usage();

    if(params)
    {
        inputParameters.Status();
    }

inputParameters.Status();






    map<string,double> Param;



    Param["model"]=model;
    Param["condition"]=condition;
    Param["phenoStrat"]=phenoStrat;
    Param["sexStrat"]=sexStrat;
    Param["parentEffect"]=parentEffect;
    Param["baseHazard"]=baseHazard;


    StatAnalysis stat(pedfile,mapfile,condpedfile,Param);
    stat.PerformAnalysis();




	return 0;

}


void vcfVersion()
{
    std::cerr << "\n\n Version  : " << VERSION<<endl;
    std::cerr << " Built    : " << DATE <<endl;
    std::cerr << " Author   : "<< USER << std::endl << std::endl;
}

void description()
{
    std::cerr <<" *******************************************************\n";
    std::cerr <<"       ModSeg - Modified Segregation Analysis \n";
    std::cerr <<"                accounting for Ascertainment Effect \n";
    std::cerr <<"                and Parent of Origin Effect" << std::endl;
    std::cerr <<" *******************************************************\n\n\n";
}

void usage()
{
    vcfVersion();
    description();
    std::cerr << "./ModSeg --ped <input pedigree file> --map <input locus file> [--out <output file>] [--params]"<< std::endl;
    std::cerr << "\t Required Parameters:\n"
              << "\t\t--ped             : Input Mendel pedigree file to read\n"
              << "\t\t--map             : Input Mendel locus file to read\n"
              << "\t\t--model           : Analysis Option to implement [Default = 0]\n"
              << "\t\t                   [Option: 0 - Case/Control Study]\n"
              << "\t\t                   [Option: 1 - Age of Onset (Survival) Study]\n\n"
              << "\t Optional Parameters (within square brackets []) :\n"
              << "\t\t--out             : Output File (Default: ModSeg.out) \n"
              << "\t\t--condition            : Conditional Analysis for Acertainment Correction (0 : No Correction [Default], 1 : Correction) \n"
              << "\t\t--sexStrat        : Stratify by Sex (0 : No Stratification [Default], 1 : Stratification)\n"
              << "\t\t--phenoStrat      : Stratify Phenotype into groups (Required if --model = 1) \n"
              << "\t\t                   [Option: 0 - No Stratificatioj]\n"
              << "\t\t                   [Option: 1 - Startify Phenotype into two equal sized groups based on boundary values in map file]\n"
              << "\t\t--condped         : Input Conditional Pedigree File to read (Required if --cond = 1) \n"
              << "\t\t--params          : Print the parameter settings\n"
              << std::endl;
}


