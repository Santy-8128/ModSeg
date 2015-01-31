


#ifndef __PEDIGREE_DATA_H
#define __PEDIGREE_DATA_H

#include <ctime>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include "Likelihood.h"
#include "Optimization.h"
#include "StringBasics.h"
#include "Pedigree.h"
#include <vector>
#include <map>
using namespace std;



template <class Ped>
class PedigreeData
{
    public:

        int DataModel;
        // 0 = for case-control study, No Age Strat, No Sex Strat
        // 1 = for survival analysis study, No Age Strat, No Sex Strat
        // 2 = for case-control study, No Age Strat, ONLY Sex Strat
        // 3 = for survival analysis study, No Age Strat, ONLY Sex Strat
        // 5 = for survival analysis study, ONLY Age Strat, No Sex Strat
        // 0 = for case-control study, No Age Strat, No Sex Strat
        // 0 = for case-control study, No Age Strat, No Sex Strat

        string MarkerName;
        string Allele1,Allele2;
        string PhenoName;
        double alleleFreq1,baseline;
        int censorIndicator;
        double phenoMin,phenoMax,phenoMid;
        vector<vector<string> > peddata;

        int NoOfPheno;
        int TotalSample;
        int NoOfPedigree;
        vector<string> PedigreeId;



        vector<Pedigree> Data;
        PedigreeData(String file,String mapfile,int model,double baselineHazard);
        PedigreeData()
        {

        };
        void ImportData();

        Ped operator[](Ped x);
        void printMapFileSummary();
        void readMapFile(String mapfile);
        Ped operator()(vector<Ped> x);
        //double FindMLE();


    //friend class Likelihood;
    friend class Node;
};




template <class Ped>
Ped PedigreeData<Ped>::operator()(vector<Ped> x)
{

    Ped lik=0.0;
    Likelihood<Ped> ThisDataLikelihood(x,alleleFreq1,DataModel,phenoMin,phenoMax,phenoMid,baseline);

    for(int i=0;i<NoOfPedigree;i++)
    {
        lik+=ThisDataLikelihood.calLogLikelihood(Data[i]);

    }

    return (0.0-lik);

}



template <class Ped>
Ped PedigreeData<Ped>::operator[](Ped x)
{

    Ped lik=0.0;
    vector<Ped> temp;
    temp.push_back(x);
    Likelihood<Ped> ThisDataLikelihood(temp,alleleFreq1,DataModel,phenoMin,phenoMax,phenoMid,baseline);

    for(int i=0;i<NoOfPedigree;i++)
    {
        lik+=ThisDataLikelihood.calLogLikelihood(Data[i]);

    }

    return (0.0-lik);

}




template <class Ped>
void PedigreeData<Ped>::ImportData()
{



    vector<vector<string> > ped;
    cout<<"\n Starting to Import Data ...\n"<<endl;

    for(int i=0;i<TotalSample;i++)
    {

                vector<string> temp;
        for(int j=1;j<(6+NoOfPheno);j++)
            temp.push_back(peddata[i][j]);
        ped.push_back(temp);
        if(i==TotalSample-1 || peddata[i+1][0].compare(peddata[i][0])!=0)
        {
           // cout<<i<<" "<<ped.size()<<" \n";
            Pedigree MyPed(ped,DataModel,Allele1,Allele2,alleleFreq1,censorIndicator,phenoMin,phenoMax);
            Data.push_back(MyPed);
            cout<<" Imported Data for Pedigree ID : "<<peddata[i][0]<<" ..."<<endl;

            PedigreeId.push_back(peddata[i][0]);
            ped.resize(0);
        }
    }





    NoOfPedigree=(int)PedigreeId.size();

    cout<<"\n Imported Data for "<<NoOfPedigree<<" pedigrees succesfully ..."<<endl;

    cout<<"\n\t  Summary of Pedigrees\n";
    cout<<"\t ----------------------\n\n";
    cout<<" No.\t Pedigree ID\t No. of Samples  \t No. of Samples\n";
    cout<<"    \t            \t [in Data]       \t [Imported]\n";
    cout<<"----\t------------\t-----------------\t-----------------\n";
    for(int i=0;i<NoOfPedigree;i++)
    {
        cout<<" "<<i+1<<"\t     "<<PedigreeId[i]<<"  \t     "<<Data[i].info.size()-1<<"  \t\t     "<<Data[i].NoSample()<<endl;
    }







}



template <class Ped>
void PedigreeData<Ped>::readMapFile(String mapfile)
{
    ifstream ifs(mapfile);


    if ( ! ifs.is_open() )
    {
        cout<< " ERROR: Cannot open Map File " << mapfile << std::endl;
        cout<< " Program Abborting ... \n\n\n\n";
        abort();
    }

    int readGeno=0;
    string line;
    int n=0;
    while(getline(ifs,line))
    {

        stringstream linestream(line);
        string data;
        getline(linestream,data,',');
        if(readGeno==2)
        {
            Allele2=data;
            readGeno=0;
        }
        else if(readGeno==1)
        {
            Allele1=data;
            getline(linestream,data,',');
            istringstream buffer(data);
            buffer>>alleleFreq1;
            readGeno=2;
        }

        else if(data.compare("GENOTYPE")==0)
        {

            getline(linestream,data,',');
            MarkerName=data;
            getline(linestream,data,',');
            int temp;
            istringstream buffer(data);
            buffer>>temp;
            if(temp!=2)
            {
                cout<< " ERROR: The Marker "<<MarkerName<<" can have only 2 alleles. \n";
                cout<< " Program Abborting ... \n\n\n\n";
                abort();
            }
            readGeno=1;
        }
        else if(data.compare("CENSOR")==0)
        {

            getline(linestream,data,',');
            istringstream buffer(data);
            buffer>>censorIndicator;
        }
        else if(data.compare("PHENOTYPE")==0)
        {

            getline(linestream,data,',');
            PhenoName=data;
            getline(linestream,data,',');
            istringstream buffer(data);
            buffer>>phenoMin;
            getline(linestream,data,',');
            istringstream buffer1(data);
            buffer1>>phenoMax;
            getline(linestream,data,',');
            istringstream buffer2(data);
            buffer2>>phenoMid;
        }
        else
        {
            cout<<" ERROR: Cannot recognize keyword \'"<<data<<"\' in Map File "<<mapfile<<endl;
            cout<<" Try --help for more details on Map Files ...\n";
            cout<<" Program Abborting ... \n\n\n\n";
            abort();

        }
    }
    ifs.close();



}


template <class Ped>
void PedigreeData<Ped>::printMapFileSummary()
{

    cout<<" Successfully Read Structure from Map File \n";
    if(DataModel==0 || DataModel==2 || DataModel==6)
        cout<<"\n\n Model Analysis Option               : Case-Control Study.\n";
    if(DataModel==1 || DataModel==3 || DataModel==5 || DataModel==7)
        cout<<" Model Analysis Option               : Age of Onset (Survival Study).\n";
    cout<<" Marker Name                         : "<<MarkerName<<endl;
    cout<<" Number of Alleles                   : 2 \n";
    cout<<" Allele 1 (Frequency)                : "<<Allele1<<" ("<< alleleFreq1 <<")\n";
    cout<<" Allele 2 (Frequency)                : "<<Allele2<<" ("<< 1-alleleFreq1 <<")\n";
    if(DataModel==1 || DataModel==3 || DataModel==5 || DataModel==7)
        {
            cout<<" Phenotype Name (Min Value,Max Value : "<<PhenoName<<" ("<<phenoMin<<","<<phenoMax<<") \n";
            cout<<" Censoring Indicator                 : "<<censorIndicator<<endl;
        }


}

template <class Ped>
PedigreeData<Ped>::PedigreeData(String pedfile,String mapfile,int model,double baselineHazard)
{
    //NoOfPheno=4;
    baseline=baselineHazard;
    DataModel=model;
    cout<<" Reading Structure from Map File "<<mapfile<<" ... "<<endl;
    readMapFile(mapfile);
    printMapFileSummary();

    if(DataModel==0 || DataModel==2 || DataModel==6)
        NoOfPheno=2;
    if(DataModel==1|| DataModel==3 || DataModel==5 || DataModel==7)
        NoOfPheno=3;



    cout<<"\n Reading Data from Pedigree File "<<pedfile<<" ..."<<endl;

    ifstream ifs(pedfile);

    if ( ! ifs.is_open() )
    {
        cout <<" ERROR: Cannot open Pedigree File " << pedfile << std::endl;
        cout<< " Program Abborting ... \n\n\n\n";
        abort();
    }

    string line;
    int n=0;
    while(getline(ifs,line))
    {

        stringstream linestream(line);
        string data;
        vector<string> temp;
        for(int i=0;i<(6+NoOfPheno);i++)
        {
            getline(linestream,data,'\t');
            //cout<<data.c_str()<<"\t";
            temp.push_back(data.c_str());
        }
        peddata.push_back(temp);
        n++;
    }
    ifs.close();
    TotalSample=(int)peddata.size();
    NoOfPedigree=0;

    cout<<" Total Number of Samples = "<<TotalSample<<endl;





}




#endif
