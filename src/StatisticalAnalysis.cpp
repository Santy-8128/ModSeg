

#include <sstream>
#include <ctime>
#include <fstream>
#include<algorithm>
#include <iostream>
#include <boost/math/distributions/normal.hpp>
#include "Likelihood.h"
#include "Pedigree.h"
#include "string.h"
# include "StatisticalAnalysis.h"
#include "Optimization.h"
#include "Numerical.h"
#include "PedigreeData.h"
#include "AscertainmentPedigreeData.h"

using namespace std;

using namespace boost::math;

double maxVal(double a,double b)
{
    return (a>=b?a:b);
}


void StatAnalysis::InvInformation()
{
    VarianceMLE.resize(InformationMatrix.size());
    singular=false;
    if((int)InformationMatrix.size()==1)
        VarianceMLE[0].push_back(1/InformationMatrix[0][0]);

    else
    {
        double det=(InformationMatrix[0][0]*InformationMatrix[1][1])-
                    (InformationMatrix[0][1]*InformationMatrix[1][0]);
       // cout<<" \n\n\n\n\n\n\n\nDET = "<<det<<endl;

        if(det<1e-4)
        {
            singular=true;
            return;
        }
        VarianceMLE[0].push_back(InformationMatrix[1][1]/det);
        VarianceMLE[0].push_back(-InformationMatrix[0][1]/det);
        VarianceMLE[1].push_back(-InformationMatrix[1][0]/det);
        VarianceMLE[1].push_back(InformationMatrix[0][0]/det);
        }



}

void StatAnalysis::printCumHazard()
{



    if(model==1 || model==3 || model==5 || model==7)
    {


        for(int i=0;i<110;i++)
            baselineHazard.push_back(base_constant);

        double sum=0.0;
        for(int i=0;i<110;i++)
        {
            sum+=baselineHazard[i];
            cummulativeBaselineHazard.push_back(sum);
        }

        vector<double> SDLambda(8);
        double SDLambdaT,SDLambdaTm,SDLambdaTf;

        cout<<"\n Estimates of Cumulative Hazard : \n";
        cout<<" ---------------------------- \n\n";

        double i=pMin;

        double SD=sqrt(VarianceMLE[0][0]);


        if(model==1)
        {
            cout<<" Phenotype \t Non-Carrier \t Carrier \t 95\% C.I. for Carrier \n";
            cout<<" ----------\t-------------\t---------\t-----------------------\n";


            while(i<=pMax)
            {

                double EstLambdaT=cummulativeBaselineHazard[i-1]*exp(MLE[0]);
                SDLambdaT=(base_constant*(i)*exp(MLE[0]))*SD;

                cout<<"  "<<i<<"\t\t"<<cummulativeBaselineHazard[i-1]<<"\t\t"<<EstLambdaT;
                cout<<"     \t ("<<maxVal(0,EstLambdaT-(critic*SDLambdaT));
                cout<<", "<<EstLambdaT+(critic*SDLambdaT)<<") \n";

                i+=10;
            }
            i=pMin;

             cout<<"\n Estimates of Cumulative Risk : \n";
                cout<<" ---------------------------- \n\n";


            cout<<" Phenotype \t Non-Carrier \t Carrier \t 95\% C.I. for Carrier \n";
            cout<<" ----------\t-------------\t---------\t-----------------------\n";


            while(i<=pMax)
            {
                double EstLambdaT=cummulativeBaselineHazard[i-1]*exp(MLE[0]);
                SDLambdaT=(base_constant*(i)*exp(MLE[0]))*SD;
                double EstRisk=1-exp(-EstLambdaT);

                cout<<"  "<<i<<"\t\t"<<1-exp(-cummulativeBaselineHazard[i-1])<<"\t"<<EstRisk;
                cout<<"     \t ("<<1-exp(-maxVal(0,EstLambdaT-(critic*SDLambdaT)));
                cout<<", "<<1-exp(-EstLambdaT-(critic*SDLambdaT))<<") \n";

                i+=10;


            }
        }
        else
        {

            double SD;
            if(model==3)
            {
                cout<<" Phenotype \t Non-Carrier \t Carrier(M) \t 95\% C.I. for Carrier(M) \t Carrier(F) \t 95\% C.I. for Carrier(F) \n";
                cout<<" ----------\t-------------\t------------\t--------------------------\t------------\t--------------------------\n";
                double SD1=sqrt(VarianceMLE[0][0]);
                double SD2=sqrt(VarianceMLE[1][1]);

                while(i<=pMax)
                {

                    double EstLambdaTm=cummulativeBaselineHazard[i-1]*exp(MLE[0]);
                    double SDLambdaTm=(base_constant*(i)*exp(MLE[0]))*SD1;
                    double EstLambdaTf=cummulativeBaselineHazard[i-1]*exp(MLE[1]);
                    double SDLambdaTf=(base_constant*(i)*exp(MLE[1]))*SD2;

                    cout<<"  "<<i<<"\t\t"<<cummulativeBaselineHazard[i-1]<<"\t\t"<<EstLambdaTm;
                    cout<<"     \t ("<<maxVal(0,EstLambdaTm-(critic*SDLambdaTm));
                    cout<<", "<<EstLambdaTm+(critic*SDLambdaTm)<<") \t\t"<<EstLambdaTf;
                    cout<<"     \t ("<<maxVal(0,EstLambdaTf-(critic*SDLambdaTf));
                    cout<<", "<<EstLambdaTf+(critic*SDLambdaTf)<<") \n";

                    i+=10;
                }
                i=pMin;

                cout<<"\n Estimates of Cumulative Risk : \n";
                cout<<" ---------------------------- \n\n";



                cout<<" Phenotype \t Non-Carrier \t Carrier(M) \t 95\% C.I. for Carrier(M) \t Carrier(F) \t 95\% C.I. for Carrier(F) \n";
                cout<<" ----------\t-------------\t------------\t--------------------------\t------------\t--------------------------\n";

                while(i<=pMax)
                {

                    double EstLambdaTm=cummulativeBaselineHazard[i-1]*exp(MLE[0]);
                    double SDLambdaTm=(base_constant*(i)*exp(MLE[0]))*SD1;
                    double EstLambdaTf=cummulativeBaselineHazard[i-1]*exp(MLE[1]);
                    double SDLambdaTf=(base_constant*(i)*exp(MLE[1]))*SD2;
                    double EstRiskm=1-exp(-EstLambdaTm);
                    double EstRiskf=1-exp(-EstLambdaTf);

                    cout<<"  "<<i<<"\t\t"<<1-exp(-cummulativeBaselineHazard[i-1])<<"\t\t"<<EstRiskm;
                    cout<<"     \t ("<<1-exp(-maxVal(0,EstLambdaTm-(critic*SDLambdaTm)));
                    cout<<", "<<1-exp(-EstLambdaTm-(critic*SDLambdaTm))<<") \t\t"<<EstRiskf;
                    cout<<"     \t ("<<1-exp(-maxVal(0,EstLambdaTf-(critic*SDLambdaTf)));
                    cout<<", "<<1-exp(-EstLambdaTf-(critic*SDLambdaTf))<<") \n";

                    i+=10;
                }

            }
            if(model==5)
            {

                  cout<<" Phenotype \t Non-Carrier \t Carrier \t 95\% C.I. for Carrier \n";
            cout<<" ----------\t-------------\t---------\t-----------------------\n";

            SD=sqrt(VarianceMLE[0][0]);
            while(i<=pMid)
            {

                double EstLambdaT=cummulativeBaselineHazard[i-1]*exp(MLE[0]);
                SDLambdaT=(base_constant*(i)*exp(MLE[0]))*SD;

                cout<<"  "<<i<<"\t\t"<<cummulativeBaselineHazard[i-1]<<"\t\t"<<EstLambdaT;
                cout<<"     \t ("<<maxVal(0,EstLambdaT-(critic*SDLambdaT));
                cout<<", "<<EstLambdaT+(critic*SDLambdaT)<<") \n";

                i+=10;
            }
            SD=sqrt(VarianceMLE[1][1]);
            while(i<=pMax)
            {

                double c1=cummulativeBaselineHazard[pMid];
                double c2=(cummulativeBaselineHazard[i-1]-cummulativeBaselineHazard[pMid]);
                double EstLambdaT=(c1*exp(MLE[0]))+((exp(MLE[1]))*c2);
                double VDLambdaT= (c1*c1*exp(2*MLE[0])*(VarianceMLE[0][0]))
                                +(c2*c2*exp(2*MLE[1])*(VarianceMLE[1][1]))
                                +(2*c1*c2*exp(MLE[0]+MLE[1])*VarianceMLE[0][1]);
                SDLambdaT=sqrt(VDLambdaT);

                cout<<"  "<<i<<"\t\t"<<cummulativeBaselineHazard[i-1]<<"\t\t"<<EstLambdaT;
                cout<<"     \t ("<<maxVal(0,EstLambdaT-(critic*SDLambdaT));
                cout<<", "<<EstLambdaT+(critic*SDLambdaT)<<") \n";

                i+=10;
            }
            i=pMin;

             cout<<"\n Estimates of Cumulative Risk : \n";
                cout<<" ---------------------------- \n\n";


            cout<<" Phenotype \t Non-Carrier \t Carrier \t 95\% C.I. for Carrier \n";
            cout<<" ----------\t-------------\t---------\t-----------------------\n";

            SD=sqrt(VarianceMLE[0][0]);
            while(i<=pMid)
            {
                double EstLambdaT=cummulativeBaselineHazard[i-1]*exp(MLE[0]);
                SDLambdaT=(base_constant*(i)*exp(MLE[0]))*SD;
                double EstRisk=1-exp(-EstLambdaT);

                cout<<"  "<<i<<"\t\t"<<1-exp(-cummulativeBaselineHazard[i-1])<<"\t\t"<<EstRisk;
                cout<<"     \t ("<<1-exp(-maxVal(0,EstLambdaT-(critic*SDLambdaT)));
                cout<<", "<<1-exp(-EstLambdaT-(critic*SDLambdaT))<<") \n";

                i+=10;


            }
            //cout<<i<<endl;
            SD=sqrt(VarianceMLE[1][1]);
            while(i<=pMax)
            {
                double c1=cummulativeBaselineHazard[pMid];
                double c2=(cummulativeBaselineHazard[i-1]-cummulativeBaselineHazard[pMid]);
                //cout<<c2<<endl;
                double EstLambdaT=(c1*exp(MLE[0]))+((exp(MLE[1]))*c2);
                double VDLambdaT= (c1*c1*exp(2*MLE[0])*(VarianceMLE[0][0]))
                                +(c2*c2*exp(2*MLE[1])*(VarianceMLE[1][1]))
                                +(2*c1*c2*exp(MLE[0]+MLE[1])*VarianceMLE[0][1]);
                SDLambdaT=sqrt(VDLambdaT);
                double EstRisk=1-exp(-EstLambdaT);


                cout<<"  "<<i<<"\t\t"<<1-exp(-cummulativeBaselineHazard[i-1])<<"\t\t"<<EstRisk;
                cout<<"     \t ("<<1-exp(-maxVal(0,EstLambdaT-(critic*SDLambdaT)));
                cout<<", "<<1-exp(-EstLambdaT-(critic*SDLambdaT))<<") \n";

//cout<<endl<<endl;
//
//                cout<<"  "<<i<<"\t\t"<<cummulativeBaselineHazard[i-1]<<"\t\t"<<EstLambdaT;
//                cout<<"     \t ("<<maxVal(0,EstLambdaT-(critic*SDLambdaT));
//                cout<<", "<<EstLambdaT+(critic*SDLambdaT)<<") \n";

                i+=10;



//
//                cout<<"  "<<i<<"\t\t"<<1-exp(-cummulativeBaselineHazard[i-1])<<"\t\t"<<EstRisk;
//                cout<<"     \t ("<<1-exp(-maxVal(0,EstLambdaT-(critic*SDLambdaT)));
//                cout<<", "<<1-exp(-EstLambdaT-(critic*SDLambdaT))<<") \n";
//
//                i+=10;


            }

            }

            if(model==7)
            {

                   cout<<" Phenotype \t Non-Carrier \t Carrier(M) \t 95\% C.I. for Carrier \n";
            cout<<" ----------\t-------------\t---------\t-----------------------\n";
SD=sqrt(VarianceMLE[0][0]);


            while(i<=pMax)
            {

                double EstLambdaT=cummulativeBaselineHazard[i-1]*exp(MLE[0]);
                SDLambdaT=(base_constant*(i)*exp(MLE[0]))*SD;

                cout<<"  "<<i<<"\t\t"<<cummulativeBaselineHazard[i-1]<<"\t\t"<<EstLambdaT;
                cout<<"     \t ("<<maxVal(0,EstLambdaT-(critic*SDLambdaT));
                cout<<", "<<EstLambdaT+(critic*SDLambdaT)<<") \n";

                i+=10;
            }
            i=pMin;

             cout<<"\n Estimates of Cumulative Risk : \n";
                cout<<" ---------------------------- \n\n";


            cout<<" Phenotype \t Non-Carrier \t Carrier(M) \t 95\% C.I. for Carrier \n";
            cout<<" ----------\t-------------\t---------\t-----------------------\n";


            while(i<=pMax)
            {
                double EstLambdaT=cummulativeBaselineHazard[i-1]*exp(MLE[0]);
                SDLambdaT=(base_constant*(i)*exp(MLE[0]))*SD;
                double EstRisk=1-exp(-EstLambdaT);

                cout<<"  "<<i<<"\t\t"<<1-exp(-cummulativeBaselineHazard[i-1])<<"\t"<<EstRisk;
                cout<<"     \t ("<<1-exp(-maxVal(0,EstLambdaT-(critic*SDLambdaT)));
                cout<<", "<<1-exp(-EstLambdaT-(critic*SDLambdaT))<<") \n";

                i+=10;

//
//
// cout<<" Phenotype \t Non-Carrier \t Carrier(M) \t 95\% C.I. for Carrier(M) \t Carrier(F) \t 95\% C.I. for Carrier(F) \n";
//                cout<<" ----------\t-------------\t------------\t--------------------------\t------------\t--------------------------\n";
//                double SD1=sqrt(VarianceMLE[0][0]);
//                double SD2=sqrt(VarianceMLE[1][1]);
//                MLE.push_
//                while(i<=pMax)
//                {
//
//                    double EstLambdaTm=cummulativeBaselineHazard[i-1]*exp(MLE[0]);
//                    double SDLambdaTm=(base_constant*(i)*exp(MLE[0]))*SD1;
//                    double EstLambdaTf=cummulativeBaselineHazard[i-1]*exp(MLE[1]);
//                    double SDLambdaTf=(base_constant*(i)*exp(MLE[1]))*SD2;
//
//                    cout<<"  "<<i<<"\t\t"<<cummulativeBaselineHazard[i-1]<<"\t\t"<<EstLambdaTm;
//                    cout<<"     \t ("<<maxVal(0,EstLambdaTm-(critic*SDLambdaTm));
//                    cout<<", "<<EstLambdaTm+(critic*SDLambdaTm)<<") \t\t"<<EstLambdaTf;
//                    cout<<"     \t ("<<maxVal(0,EstLambdaTf-(critic*SDLambdaTf));
//                    cout<<", "<<EstLambdaTf+(critic*SDLambdaTf)<<") \n";
//
//                    i+=10;
//                }
//                i=pMin;
//
//                cout<<"\n Estimates of Cumulative Risk : \n";
//                cout<<" ---------------------------- \n\n";
//
//
//
//                cout<<" Phenotype \t Non-Carrier \t Carrier(M) \t 95\% C.I. for Carrier(M) \t Carrier(F) \t 95\% C.I. for Carrier(F) \n";
//                cout<<" ----------\t-------------\t------------\t--------------------------\t------------\t--------------------------\n";
//
//                while(i<=pMax)
//                {
//
//                    double EstLambdaTm=cummulativeBaselineHazard[i-1]*exp(MLE[0]);
//                    double SDLambdaTm=(base_constant*(i)*exp(MLE[0]))*SD1;
//                    double EstLambdaTf=cummulativeBaselineHazard[i-1]*exp(MLE[1]);
//                    double SDLambdaTf=(base_constant*(i)*exp(MLE[1]))*SD2;
//                    double EstRiskm=1-exp(-EstLambdaTm);
//                    double EstRiskf=1-exp(-EstLambdaTf);
//
//                    cout<<"  "<<i<<"\t\t"<<1-exp(-cummulativeBaselineHazard[i-1])<<"\t\t"<<EstRiskm;
//                    cout<<"     \t ("<<1-exp(-maxVal(0,EstLambdaTm-(critic*SDLambdaTm)));
//                    cout<<", "<<1-exp(-EstLambdaTm-(critic*SDLambdaTm))<<" \t\t"<<EstRiskf;
//                    cout<<"     \t ("<<1-exp(-maxVal(0,EstLambdaTf-(critic*SDLambdaTf)));
//                    cout<<", "<<1-exp(-EstLambdaTf-(critic*SDLambdaTf))<<") \n";
//
//                    i+=10;
//                }
            }



        }


    }




}
}




void StatAnalysis::printCumRisk()
{





}


void StatAnalysis::resultPrint()
{

    normal norm(0,1);

    cout<<" Maximized Log-Likelihood is : "<<-MaxLogLikelihood<<endl<<endl;

    cout<<" Maximum Likelihood Estimates  : \n";
    cout<<" ---------------------------- \n\n";


    if(model==1 || model==3 || model==5 || model==7)
        cout<<"\n Estimate of Log of Hazard Ratio\n\n";
    if(model==2 || model==4 || model==6)
        cout<<"\n Estimate of Effect\n";

    if(singular)
    {

        cout<<" Parameter \t\t\t Estimate  \n";
        cout<<" ----------\t\t\t----------\n";
        for(int i=0;i<(int)MLE.size();i++)
            cout<<ParamName[i]<<"\t\t "<<MLE[i]<<endl;

        cout<<" WARNING !!! Estimated Variance-Covariance Matrix is singular !!! P-values could NOT be obtained !!! \n\n";
        abort();
    }

    for(int i=0;i<(int)MLE.size();i++)
    {

        mlePvalue.push_back(2*cdf(complement(norm,abs(MLE[i]/sqrt(VarianceMLE[i][i])))));
    }


    cout<<" Parameter \t\t\t Estimate \t S.E. \t\t 95\% C.I. \t\t P-Value \n";
    cout<<" ----------\t\t\t----------\t------\t\t-----------\t\t---------\n";

    critic=1.96;


    for(int i=0;i<(int)MLE.size();i++)
    {
        cout<<ParamName[i]<<"\t\t "<<MLE[i]<<"   \t"<<sqrt(VarianceMLE[i][i])<<"     \t";
        cout<<"("<<MLE[i]-(critic*sqrt(VarianceMLE[i][i]))<<", "<<MLE[i]+(critic*sqrt(VarianceMLE[i][i]))<<")   \t";
        cout<<mlePvalue[i]<<endl;
    }

    cout<<"\n Estimated Fisher's Information Matrix : \n";
    cout<<" ------------------------------------- \n\n";


    for(int i=0;i<(int)MLE.size();i++)
    {
        cout<<"\t";
        for(int j=0;j<(int)MLE.size();j++)
        {

            cout<<InformationMatrix[i][j]<<"\t\t";

        }
        cout<<endl;

    }


    printCumHazard();

    cout<<endl;


}



StatAnalysis::StatAnalysis(String pedfile,String mapfile,String condPedfile,map<string,double> &Param)
{

    pedFile=pedfile;
    condPedFile=condPedfile;
    mapFile=mapfile;
    pMid=40;
    base_constant=0.0001;


    map<string,double>::iterator it;
    for(it=Param.begin();it!=Param.end();++it)
    {

        model=0;
        if(Param.find("model")!=Param.end())
            model=Param.find("model")->second;


        if(model!=0 && model!=1)
        {
            cout <<" ERROR: Invalid --model Parameter = "<<Param.find("model")->second <<std::endl;
            cout <<" See Usage for more details ...\n";
            cout<< " Program Abborting ... \n\n\n\n";
            abort();
        }


        condition=0;
        if(Param.find("condition")!=Param.end())
            condition=Param.find("condition")->second;
        if(condition!=0 && condition!=1)
        {
            cout <<" ERROR: Invalid --condition Parameter = "<<Param.find("condition")->second <<std::endl;
            cout <<" See Usage for more details ...\n";
            cout<< " Program Abborting ... \n\n\n\n";
            abort();
        }

        phenoStrat=0;
        if(Param.find("phenoStrat")!=Param.end())
            phenoStrat=Param.find("phenoStrat")->second;
        if(phenoStrat!=0 && phenoStrat!=1)
        {
            cout <<" ERROR: Invalid --phenoStrat Parameter = "<<Param.find("phenoStrat")->second <<std::endl;
            cout <<" See Usage for more details ...\n";
            cout<< " Program Abborting ... \n\n\n\n";
            abort();
        }

        sexStrat=0;
        if(Param.find("sexStrat")!=Param.end())
            sexStrat=Param.find("sexStrat")->second;
        if(sexStrat!=0 && sexStrat!=1)
        {
            cout <<" ERROR: Invalid --sexStrat Parameter = "<<Param.find("sexStrat")->second <<std::endl;
            cout <<" See Usage for more details ...\n";
            cout<< " Program Abborting ... \n\n\n\n";
            abort();
        }



        parentEffect=0;
        if(Param.find("parentEffect")!=Param.end())
            parentEffect=Param.find("parentEffect")->second;
        if(parentEffect!=0 && parentEffect!=1)
        {
            cout <<" ERROR: Invalid --parentEffect Parameter = "<<Param.find("parentEffect")->second <<std::endl;
            cout <<" See Usage for more details ...\n";
            cout<< " Program Abborting ... \n\n\n\n";
            abort();
        }


//        if(Param.find("baseHazard")!=Param.end())
//            base_constant=Param.find("baseHazard")->second;


        if(phenoStrat==1 && sexStrat==1)
        {
            cout <<" ERROR: Cannot stratify by both Phenotype and Sex. "<<std::endl;
            cout<< " Only one of --phenoStrat and --sexStrat can be 1 (too many parameters) !!!"<<endl;
            cout<< " Program Abborting ... \n\n\n\n";
            abort();
        }

        if(parentEffect==1 && sexStrat==1)
        {
            cout <<" ERROR: Cannot have Parent of Origin Effect as well as stratify by Sex. "<<std::endl;
            cout<< " Only one of --parentEffect and --sexStrat can be 1 (too many parameters) !!!"<<endl;
            cout<< " Program Abborting ... \n\n\n\n";
            abort();
        }

        if(parentEffect==1 && phenoStrat==1)
        {
            cout <<" ERROR: Cannot have Parent of Origin Effect as well as stratify by Phenotype. "<<std::endl;
            cout<< " Only one of --parentEffect and --phenoStrat can be 1 (too many parameters) !!!"<<endl;
            cout<< " Program Abborting ... \n\n\n\n";
            abort();
        }


        if(phenoStrat==1 && model==0)
        {
            cout <<" ERROR: Cannot stratify Phenotype for Case-Control Study. "<<std::endl;
            cout<< " Phenotype should be continuous for stratification."<<endl;
            cout<< " Program Abborting ... \n\n\n\n";
            abort();
        }

        if(condition==1 && condPedfile=="")
        {
            cout <<" ERROR: No Conditional Pedigree File provided for correction of ascertainment bias. "<<std::endl;
            cout<< " See Usage for more details ..."<<endl;
            cout<< " Program Abborting ... \n\n\n\n";
            abort();

        }


    }
}

void StatAnalysis::PerformAnalysis()
{

    int tempModel=model;
    ParamName.push_back(" First Allele       ");
    if(tempModel==1)
    {
        if(sexStrat==0 && phenoStrat==0 && parentEffect==0)
            model=1;
        else if(sexStrat==1)
        {
            model=3;
            ParamName.push_back(" Female             ");
        }
        else if(phenoStrat==1)
        {
            model=5;
            ParamName.push_back(" Second Phenotype Group");
        }
        else
        {
            model=7;
            //ParamName.push_back(" Parent Effect (Mother)");
        }
    }
    else if(tempModel==0)
    {
        if(sexStrat==0 && parentEffect==0)
            model=0;
        else if(sexStrat==1)
        {
            model=2;
            ParamName.push_back(" Female             ");
        }
        else
        {
            model=6;
            ParamName.push_back(" Parent of Origin Effect (Mother)");
        }
    }


    if(condition==0)
    {


        PedigreeData<double> ThisData(pedFile,mapFile,model,base_constant);
        ThisData.ImportData();
        cout<<"\n\n Calculating Maximum Likelihood Estimates ...\n\n";

         if(parentEffect==1)
        {
            double x[1];
            x[0]=0.0;
            Optimization<PedigreeData<double> > ThisDataMLE(x,1);
//            ThisDataMLE.amoeba(ThisData,1e-20);
//            MLE=ThisDataMLE.xmin();
               vector<double> temp(1);
            temp[0]=0.0;
            //temp[1]=0.0;

            //cout<<" VALUE = -"<<ThisData(temp)<<endl;

             for(double ll=0;ll<=10;ll=ll+0.1)
            {
                temp[0]=ll;
                cout<<" VALUE = "<<ll<<" \t"<<ThisData(temp)<<endl;


            }




            ThisDataMLE.amoeba(ThisData,1e-20);
            MLE=ThisDataMLE.xmin();

            cout<<" FINAL VALUE CHECK = "<<MLE[0]<<" "<<ThisData(MLE)<<endl;




            MaxLogLikelihood=ThisDataMLE.ymin();

        }
        else if(phenoStrat==0 && sexStrat==0)
        {

            double x[1];
            x[0]=0.0;
            Optimization<PedigreeData<double> > ThisDataMLE(x,1);
//            ThisDataMLE.amoeba(ThisData,1e-20);
//            MLE=ThisDataMLE.xmin();

                vector<double> temp(1);
            temp[0]=0.0;
            //temp[1]=0.0;

            //cout<<" VALUE = -"<<ThisData(temp)<<endl;

             for(double ll=0;ll<=10;ll=ll+0.1)
            {
                temp[0]=ll;
                cout<<" VALUE = "<<ll<<" \t"<<ThisData(temp)<<endl;


            }




            ThisDataMLE.amoeba(ThisData,1e-20);
            MLE=ThisDataMLE.xmin();

            cout<<" FINAL VALUE CHECK = "<<MLE[0]<<" "<<ThisData(MLE)<<endl;




            MaxLogLikelihood=ThisDataMLE.ymin();
        }
        else
        {
            double x[2];
            x[0]=0.0;
            x[1]=0.0;
            Optimization<PedigreeData<double> > ThisDataMLE(x,2);
            ThisDataMLE.amoeba(ThisData,1e-20);
            MLE=ThisDataMLE.xmin();
            MaxLogLikelihood=ThisDataMLE.ymin();

        }

        cout<<"\n Calculating Fisher's Information Matrix ...\n\n";

        SecondDifferentiation<PedigreeData<double> > FisherInfo;
        InformationMatrix=FisherInfo.doubleDifferential(ThisData,MLE,1e-3);
        InvInformation();
         if(model==1 || model==3 || model ==5 || model ==7)
                {
                    pMin=ThisData.phenoMin;
                    pMax=ThisData.phenoMax;
                    pMid=ThisData.phenoMid;
                }


    }
    else
    {

        AscertainmentPedigreeData ThisData(pedFile,condPedFile,mapFile,model,base_constant);
        ThisData.ImportData();
        cout<<"\n\n Calculating Maximum Likelihood Estimates ...\n\n";

        if(parentEffect==1)
        {
            double x[1];
            x[0]=0.0;
            Optimization<AscertainmentPedigreeData> ThisDataMLE(x,1);
//            ThisDataMLE.amoeba(ThisData,1e-20);
//            MLE=ThisDataMLE.xmin();

            vector<double> temp(1);
            temp[0]=0.0;
            //temp[1]=0.0;

            //cout<<" VALUE = -"<<ThisData(temp)<<endl;

             for(double ll=0;ll<=10;ll=ll+0.1)
            {
                temp[0]=ll;
                cout<<" VALUE = "<<ll<<" \t"<<ThisData(temp)<<endl;


            }




            ThisDataMLE.amoeba(ThisData,1e-20);
            MLE=ThisDataMLE.xmin();

            cout<<" FINAL VALUE CHECK = "<<MLE[0]<<" "<<ThisData(MLE)<<endl;


            MaxLogLikelihood=ThisDataMLE.ymin();

        }
        else if(phenoStrat==0 && sexStrat==0)
        {
            double x[1];
            x[0]=0.0;
            Optimization<AscertainmentPedigreeData> ThisDataMLE(x,1);
//            ThisDataMLE.amoeba(ThisData,1e-20);
//            MLE=ThisDataMLE.xmin();

               vector<double> temp(1);
            temp[0]=0.0;
            //temp[1]=0.0;

            //cout<<" VALUE = -"<<ThisData(temp)<<endl;

             for(double ll=0;ll<=10;ll=ll+0.1)
            {
                temp[0]=ll;
                cout<<" VALUE = "<<ll<<" \t"<<ThisData(temp)<<endl;


            }




            ThisDataMLE.amoeba(ThisData,1e-20);
            MLE=ThisDataMLE.xmin();

            cout<<" FINAL VALUE CHECK = "<<MLE[0]<<" "<<ThisData(MLE)<<endl;




            MaxLogLikelihood=ThisDataMLE.ymin();
        }
        else
        {
            double x[2];
            x[0]=0.0;
            x[1]=0.0;
            Optimization<AscertainmentPedigreeData> ThisDataMLE(x,2);
            ThisDataMLE.amoeba(ThisData,1e-20);



MLE=ThisDataMLE.xmin();
            MaxLogLikelihood=ThisDataMLE.ymin();

        }
        cout<<"\n Calculating Fisher's Information Matrix ...\n\n";

        SecondDifferentiation<AscertainmentPedigreeData> FisherInfo;
        InformationMatrix=FisherInfo.doubleDifferential(ThisData,MLE,1e-2);
        InvInformation();
        if(model==1 || model==3 || model ==5 || model ==7)
        {
            pMin=ThisData.phenoMin;
            pMax=ThisData.phenoMax;
        }
    }



    cout<<" Printing out Results ... \n\n";

    resultPrint();






}


