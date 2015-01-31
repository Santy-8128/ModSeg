#include<algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <ctime>
#include "Likelihood.h"
#include "Optimization.h"
#include "AscertainmentPedigreeData.h"
//#include "string.h"


using namespace std;


 double AscertainmentPedigreeData::operator()(vector<double> x)
{

    double lik=0.0;
    Likelihood<double> ThisDataLikelihood(x,alleleFreq1,DataModel,phenoMin,phenoMax,phenoMid,baseline);

    for(int i=0;i<NoOfPedigree1;i++)
    {
        lik+=ThisDataLikelihood.calLogLikelihood(Data[i]);
    }


    for(int i=0;i<NoOfPedigree2;i++)
    {
        lik-=ThisDataLikelihood.calLogLikelihood(CondData[i]);
    }


    return (0.0-lik);

}



 double AscertainmentPedigreeData::operator[](double x)
{

    double lik=0.0;
    vector<double> temp;
    temp.push_back(x);
    Likelihood<double> ThisDataLikelihood(temp,alleleFreq1,DataModel,phenoMin,phenoMax,phenoMid,baseline);

    for(int i=0;i<NoOfPedigree1;i++)
    {
        lik+=ThisDataLikelihood.calLogLikelihood(Data[i]);
    }


    for(int i=0;i<NoOfPedigree2;i++)
    {
        lik-=ThisDataLikelihood.calLogLikelihood(CondData[i]);
    }


    return (0.0-lik);

}













void AscertainmentPedigreeData::ImportData()
{
    vector<vector<string> > ped;
    cout<<"\n Starting to Import Data for Main Pedigree...\n"<<endl;

    for(int i=0;i<TotalSample1;i++)
    {

        vector<string> temp;
        for(int j=1;j<(6+NoOfPheno);j++)
            temp.push_back(peddata[i][j]);
        ped.push_back(temp);
        if(i==TotalSample1-1 || peddata[i+1][0].compare(peddata[i][0])!=0)
        {
            Pedigree MyPed(ped,DataModel,Allele1,Allele2,alleleFreq1,censorIndicator,phenoMin,phenoMax);
            Data.push_back(MyPed);
            PedigreeId1.push_back(peddata[i][0]);
            ped.resize(0);
        }
    }

    NoOfPedigree1=(int)PedigreeId1.size();

    cout<<"\n Imported Data for "<<NoOfPedigree1<<" pedigrees succesfully ..."<<endl;

    cout<<"\n\t  Summary of Main Pedigrees\n";
    cout<<"\t ----------------------------\n\n";
    cout<<" No.\t Pedigree ID\t No. of Samples  \t No. of Samples\n";
    cout<<"    \t            \t [in Data]       \t [Imported]\n";
    cout<<"----\t------------\t-----------------\t-----------------\n";
    for(int i=0;i<NoOfPedigree1;i++)
    {
        cout<<" "<<i+1<<"\t     "<<PedigreeId1[i]<<"  \t     "<<Data[i].info.size()-1<<"  \t\t     "<<Data[i].NoSample()<<endl;
    }




    cout<<"\n Starting to Import Data for Conditional Pedigree...\n"<<endl;

    for(int i=0;i<TotalSample2;i++)
    {

        vector<string> temp;
        for(int j=1;j<(6+NoOfPheno);j++)
            temp.push_back(Condpeddata[i][j]);
        ped.push_back(temp);
        if(i==TotalSample1-1 || Condpeddata[i+1][0].compare(Condpeddata[i][0])!=0)
        {
            Pedigree MyPed(ped,DataModel,Allele1,Allele2,alleleFreq1,censorIndicator,phenoMin,phenoMax);
            CondData.push_back(MyPed);
            PedigreeId2.push_back(Condpeddata[i][0]);
            ped.resize(0);
        }
    }

    NoOfPedigree2=(int)PedigreeId2.size();

    cout<<"\n Imported Data for "<<NoOfPedigree2<<" pedigrees succesfully ..."<<endl;

    cout<<"\n\t  Summary of Conditional Pedigrees\n";
    cout<<"\t -----------------------------------\n\n";
    cout<<" No.\t Pedigree ID\t No. of Samples  \t No. of Samples\n";
    cout<<"    \t            \t [in Data]       \t [Imported]\n";
    cout<<"----\t------------\t-----------------\t-----------------\n";
    for(int i=0;i<NoOfPedigree2;i++)
    {
        cout<<" "<<i+1<<"\t     "<<PedigreeId2[i]<<"  \t     "<<CondData[i].info.size()-1<<"  \t\t     "<<CondData[i].NoSample()<<endl;
    }




}



void AscertainmentPedigreeData::readMapFile(char* mapfile)
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


void AscertainmentPedigreeData::printMapFileSummary()
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



void AscertainmentPedigreeData::readPedigreeFile(char* file,int cond)
{

    ifstream ifs(file);

    if ( ! ifs.is_open() )
    {
        cout <<" ERROR: Cannot open Pedigree File " << file << std::endl;
        cout<< " Program Abborting ... \n\n\n\n";
        abort();
    }
    string line;
    while(getline(ifs,line))
    {
        stringstream linestream(line);
        string data;
        vector<string> temp;
        for(int i=0;i<(6+NoOfPheno);i++)
        {
            getline(linestream,data,'\t');
            temp.push_back(data.c_str());
        }

        if(cond==0)
            peddata.push_back(temp);
        if(cond==1)
            Condpeddata.push_back(temp);

    }
    ifs.close();

}




AscertainmentPedigreeData::AscertainmentPedigreeData(char* pedfile,char* condpedfile,char* mapfile,int model,double baselineHazard)
{

baseline=baselineHazard;
    DataModel=model;
    cout<<" Reading Structure from Map File "<<mapfile<<" ... "<<endl;
    readMapFile(mapfile);
    printMapFileSummary();

   if(DataModel==0 || DataModel==2 || DataModel==4 )
        NoOfPheno=2;
    if(DataModel==1|| DataModel==3 || DataModel==5 || DataModel==7)
        NoOfPheno=3;




    cout<<"\n Reading Data from Main Pedigree File "<<pedfile<<" ..."<<endl;

    readPedigreeFile(pedfile,0);
    TotalSample1=(int)peddata.size();
    cout<<" Total Number of Samples in Main Pedigree File = "<<TotalSample1<<endl;
    NoOfPedigree1=0;


    cout<<"\n Reading Data from Conditional Pedigree File "<<pedfile<<" ..."<<endl;

    readPedigreeFile(condpedfile,1);
    TotalSample2=(int)Condpeddata.size();
    cout<<" Total Number of Samples in Main Pedigree File = "<<TotalSample2<<endl;
    NoOfPedigree2=0;


}




