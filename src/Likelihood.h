
#ifndef __LIKELIHOOD_H
#define __LIKELIHOOD_H


#include <string>
#include <vector>
#include <map>
#include<algorithm>
#include <iostream>
#include <cmath>
#include "Pedigree.h"
#include "string.h"
//#include "PedigreeData.h"

#define ZEPS 1e-310

using namespace std;




template <class Lik>
class Likelihood

{
    protected:

        vector<double> logLikePedigree;

    public:


        Lik  lambda,phenoStrat,parentEffect,sexStrat,alleleFreq;
        vector<double> baselineHazard;
        vector<double> cummulativeBaselineHazard;


        map<Node*,map<string, Lik > > LikelihoodValue;
        Lik p;
        vector<string> geno;
        int DataModel;
        double midPheno;
        Likelihood(vector<Lik> &param,Lik &AF,int &model,double phenoMin,double phenoMax,double Mp,double baseline);




        Lik PhenoGivenGeno(Phenotype &pheno,string geno,string sex,string fatherGeno,string motherGeno);
        double ChildGivenParent(string genoChild,string genoSpouse,string genoParent);
        Lik RecursionLikelihood(Node* current,Node *from,string geno,string direction);
        Lik ProbGeno(string genotype);
        Lik calLogLikelihood(Pedigree &ped);

};


template <class Lik>
Likelihood<Lik>::Likelihood(vector<Lik> &param,Lik &AF,int &model,double phenoMin,double phenoMax,double Mp,double baseline)
{

    DataModel=model;
    if(DataModel==1)
    {
        lambda=param[0];
    }

    if(DataModel==3)
    {
        lambda=param[0];
        sexStrat=param[1];
    }


    if(DataModel==5)
    {
        lambda=param[0];
        phenoStrat=param[1];
    }

    if(DataModel==7)
    {
        lambda=param[0];
        parentEffect=0.0;
    }


    midPheno=Mp;


    p=AF;
    string s0="0";
    string s1="1";
    string s2="2";

    geno.push_back(s0);
    geno.push_back(s1);
    geno.push_back(s2);

    double base_constant=baseline;
    for(int i=0;i<110;i++)
        baselineHazard.push_back(base_constant);

    double sum=0.0;
    for(int i=0;i<110;i++)
    {
        sum+=baselineHazard[i];
        cummulativeBaselineHazard.push_back(sum);
    }


}



template <class Lik>
Lik Likelihood<Lik>::ProbGeno(string genotype)
{
    if(genotype.compare(geno[0])==0)
        return (p*p);

    if(genotype.compare(geno[1])==0)
        return (2*p*(1-p));

    if(genotype.compare(geno[2])==0)
        return ((1-p)*(1-p));

}


template <class Lik>
Lik Likelihood<Lik>::PhenoGivenGeno(Phenotype &pheno,string geno,string sex,string fatherGeno,string motherGeno)
{
    //cout<<DataModel<<endl;
    if(pheno.ageOfOnset==-1)
            {
                return 1.0;
            }


    if(DataModel==1)
    {

        if(geno.compare("2")==0)
        {
            if(pheno.censor==false)
            {
                return baselineHazard[pheno.ageOfOnset]*exp(-cummulativeBaselineHazard[pheno.ageOfOnset]);
            }
            if(pheno.censor==true)
            {
                return exp(-cummulativeBaselineHazard[pheno.ageOfOnset]);
            }
        }

        else
        {
            if(pheno.censor==false)
            {
                return exp(lambda)*baselineHazard[pheno.ageOfOnset]*pow(exp(-cummulativeBaselineHazard[pheno.ageOfOnset]),exp(lambda));
            }
            if(pheno.censor==true)
            {
                return pow(exp(-cummulativeBaselineHazard[pheno.ageOfOnset]),exp(lambda));
            }
        }
    }



    else if(DataModel==3)
    {

        if(geno.compare("2")==0)
        {
            if(pheno.censor==false)
            {
                return baselineHazard[pheno.ageOfOnset]*exp(-cummulativeBaselineHazard[pheno.ageOfOnset]);
            }
            if(pheno.censor==true)
            {
                return exp(-cummulativeBaselineHazard[pheno.ageOfOnset]);
            }
        }
        else
        {
            if(sex.compare("M")==0)
            {
                if(pheno.censor==false)
                {
                    return exp(lambda)*baselineHazard[pheno.ageOfOnset]*pow(exp(-cummulativeBaselineHazard[pheno.ageOfOnset]),exp(lambda));
                }
                if(pheno.censor==true)
                {
                    return pow(exp(-cummulativeBaselineHazard[pheno.ageOfOnset]),exp(lambda));
                }
            }

            if(sex.compare("F")==0)
            {
                if(pheno.censor==false)
                {
                    return exp(sexStrat)*baselineHazard[pheno.ageOfOnset]*pow(exp(-cummulativeBaselineHazard[pheno.ageOfOnset]),exp(sexStrat));
                }
                if(pheno.censor==true)
                {
                    return pow(exp(-cummulativeBaselineHazard[pheno.ageOfOnset]),exp(sexStrat));
                }
            }

        }
    }




    else if(DataModel==5)
    {
        //
        if(geno.compare("2")==0)
        {
            if(pheno.censor==false)
            {
                return baselineHazard[pheno.ageOfOnset]*exp(-cummulativeBaselineHazard[pheno.ageOfOnset]);
            }
            if(pheno.censor==true)
            {
                return exp(-cummulativeBaselineHazard[pheno.ageOfOnset]);
            }
        }
        else
        {
            if(pheno.ageOfOnset<midPheno)
            {

                //cout<<midPheno<<"\t"<<pheno.ageOfOnset<<endl;
                if(pheno.censor==false)
                {
                    return exp(lambda)*baselineHazard[pheno.ageOfOnset]*exp(-((cummulativeBaselineHazard[pheno.ageOfOnset])*(exp(lambda))));
                }
                if(pheno.censor==true)
                {
                    return exp(-((cummulativeBaselineHazard[pheno.ageOfOnset])*(exp(lambda))));
                }
            }

            if(pheno.ageOfOnset>=midPheno)
            {

                double LL= ((cummulativeBaselineHazard[midPheno])*(exp(lambda)))
                            + ((exp(phenoStrat))*(cummulativeBaselineHazard[pheno.ageOfOnset])-(cummulativeBaselineHazard[midPheno]));
                double ll=exp(phenoStrat)*baselineHazard[pheno.ageOfOnset];


                if(pheno.censor==false)
                {
                    return ll*exp(-LL);
                }
                if(pheno.censor==true)
                {
                    return exp(-LL);
                }
            }

        }
    }



    else if(DataModel==7)
    {


        double fatherNoCensor=exp(lambda)*baselineHazard[pheno.ageOfOnset]*exp(-cummulativeBaselineHazard[pheno.ageOfOnset]*exp(lambda));;
        double fatherCensor=exp(-cummulativeBaselineHazard[pheno.ageOfOnset]*exp(lambda));

//        double fatherNoCensor=exp(parentEffect)*baselineHazard[pheno.ageOfOnset]*exp(-cummulativeBaselineHazard[pheno.ageOfOnset]*exp(parentEffect));
//        double fatherCensor=exp(-cummulativeBaselineHazard[pheno.ageOfOnset]*exp(parentEffect));
        double motherNoCensor=exp(parentEffect)*baselineHazard[pheno.ageOfOnset]*exp(-cummulativeBaselineHazard[pheno.ageOfOnset]*exp(parentEffect));
        double motherCensor=exp(-cummulativeBaselineHazard[pheno.ageOfOnset]*exp(parentEffect));


        double bothNoCensor=1-((1-fatherNoCensor)*(1-motherNoCensor));
        double bothCensor=1-((1-fatherCensor)*(1-motherCensor));


        //exp(lambda)*baselineHazard[pheno.ageOfOnset]*pow(exp(-cummulativeBaselineHazard[pheno.ageOfOnset]),exp(lambda));

        if(geno.compare("2")==0)
        {
            if(pheno.censor==false)
                return baselineHazard[pheno.ageOfOnset]*exp(-cummulativeBaselineHazard[pheno.ageOfOnset]);
            if(pheno.censor==true)
                return exp(-cummulativeBaselineHazard[pheno.ageOfOnset]);
        }
        else if(geno.compare("1")==0)
        {


            if(fatherGeno.compare("2")==0)
            {
                if(pheno.censor==false)
                    return motherNoCensor;
                if(pheno.censor==true)
                    return motherCensor;

            }
            else if(motherGeno.compare("2")==0)
            {
                if(pheno.censor==false)
                    return fatherNoCensor;
                if(pheno.censor==true)
                    return fatherCensor;

            }
            else
            {

                if(fatherGeno.compare("1")==0)
                {
                    if(motherGeno.compare("0")==0)
                    {
                        if(pheno.censor==false)
                            return motherNoCensor;
                        if(pheno.censor==true)
                            return motherCensor;
                    }

                    else
                    {
                        if(pheno.censor==false)
                            return (motherNoCensor+fatherNoCensor)/2;
                        if(pheno.censor==true)
                            return (motherCensor+fatherCensor)/2;

                    }
                }
                else
                {
                    if(pheno.censor==false)
                        return fatherNoCensor;
                    if(pheno.censor==true)
                        return fatherCensor;
                }




            }

        }

        else
        {
            if(fatherGeno.compare("0")==0 && motherGeno.compare("0")==0)
            {

            }
            else if(fatherGeno.compare("0")==0 && motherGeno.compare("1")==0)
            {

            }
            else if(fatherGeno.compare("1")==0 && motherGeno.compare("0")==0)
            {

            }
            else
            {
                if(pheno.censor==false)
                    return bothNoCensor;
                if(pheno.censor==true)
                    return bothCensor;
            }


        }
    }






    if(geno.compare("0")==0)
    {
        return 1.0;
    }
    if(geno.compare("1")==0)
    {
        return 1.0;
    }
    if(geno.compare("2")==0)
    {
        return 1.0;
    }

}

template <class Lik>
double Likelihood<Lik>::ChildGivenParent(string genoChild,string genoSpouse,string genoParent)
{


    vector<vector<vector<double> > > genoChilddata;
    genoChilddata.resize(3);

    genoChilddata[0].resize(3);
    genoChilddata[1].resize(3);
    genoChilddata[2].resize(3);

    genoChilddata[0][0].resize(3,0.0);
    genoChilddata[0][1].resize(3,0.0);
    genoChilddata[0][2].resize(3,0.0);
    genoChilddata[1][0].resize(3,0.0);
    genoChilddata[1][1].resize(3,0.0);
    genoChilddata[1][2].resize(3,0.0);
    genoChilddata[2][0].resize(3,0.0);
    genoChilddata[2][1].resize(3,0.0);
    genoChilddata[2][2].resize(3,0.0);



 //cout<<"I KNOWWWWWWW RIGHTTTTTT\n";
    genoChilddata[0][0][0]=1.0;

    genoChilddata[0][0][0]=1.0;

    genoChilddata[0][1][0]=0.5;
    genoChilddata[0][1][1]=0.5;

    genoChilddata[0][2][1]=1.0;


    genoChilddata[1][0][0]=0.5;
    genoChilddata[1][0][1]=0.5;

    genoChilddata[1][1][0]=0.25;
    genoChilddata[1][1][1]=0.5;
    genoChilddata[1][1][2]=0.25;

    genoChilddata[1][2][1]=0.5;
    genoChilddata[1][2][2]=0.5;

    genoChilddata[2][0][1]=1.0;

    genoChilddata[2][1][1]=0.5;
    genoChilddata[2][1][2]=0.5;

    genoChilddata[2][2][2]=1.0;



    for(int i=0;i<(int)geno.size();i++)
    {
        for(int j=0;j<(int)geno.size();j++)
        {
            for(int k=0;k<(int)geno.size();k++)
            {
                if(genoChild.compare(geno[k])==0 && (genoSpouse.compare(geno[j])==0 && genoParent.compare(geno[i])==0))
                    return genoChilddata[i][j][k];
            }
        }
    }







}



template <class Lik>
Lik Likelihood<Lik>::RecursionLikelihood(Node* current,Node *from,string genotype,string direction)
{

    if(LikelihoodValue.count(current)>0)
    {
        if(LikelihoodValue[current].count(genotype)>0)
            return LikelihoodValue[current][genotype];
    }

    if(direction.compare("child")==0)
    {
        double allSpouseLik=1.0;

        for(int i=0;i<(int)current->spouse.size();i++)
        {




            Node* thisSpouse=current->spouse[i];
            double spouseLik=0;

            vector<string> SpouseGeno;
            if(thisSpouse->geno.compare("")==0)
                SpouseGeno=geno;
            else
                SpouseGeno.push_back(thisSpouse->geno);



            for(int l=0;l<(int)SpouseGeno.size();l++)
            {

//                string fatherGenotype=(current->sex.compare("M")==0?current->geno:thisSpouse->geno);
//                string motherGenotype=(current->sex.compare("M")==0?thisSpouse->geno:current->geno);

                double childLik=1.0;
                for(int j=0;j<(int)current->child[thisSpouse].size();j++)
                {

                    Node* thisChild=current->child[thisSpouse][j];

                    double childGenoSum=0.0;
                    //cout<<"WELLLL\n";
                    vector<string> ChildGeno;
                    if(thisChild->geno.compare("")==0)
                        ChildGeno=geno;
                    else
                        ChildGeno.push_back(thisChild->geno);



                    for(int m=0;m<(int)ChildGeno.size();m++)
                    {
                        childGenoSum+=(PhenoGivenGeno(thisChild->pheno,ChildGeno[m],thisChild->sex,thisChild->father->geno,thisChild->mother->geno)*ChildGivenParent(ChildGeno[m],SpouseGeno[l],genotype)*RecursionLikelihood(thisChild,current,ChildGeno[m],"child"));

                    }
                    childLik*=childGenoSum;
                }
                string fatherGenotype="";
                string motherGenotype="";
                if(!thisSpouse->founder)
                {
                    fatherGenotype=thisSpouse->father->geno;
                    motherGenotype=thisSpouse->mother->geno;
                }
                spouseLik+=(childLik*PhenoGivenGeno(thisSpouse->pheno,SpouseGeno[l],thisSpouse->sex,fatherGenotype,motherGenotype)*RecursionLikelihood(thisSpouse,current,SpouseGeno[l],"spouse"));


            }
            allSpouseLik*=spouseLik;

        }
        //cout<<current->ID<<"\t"<<genotype<<"\t"<<allSpouseLik<<endl;

        LikelihoodValue[current][genotype]=allSpouseLik;

        return allSpouseLik;
    }






if(direction.compare("parent")==0)
    {
        double finalLikelihood=0.0;

        vector<string> currentGeno;
        if(current->geno.compare("")==0)
            currentGeno=geno;
        else
            currentGeno.push_back(current->geno);



        for(int k=0;k<(int)currentGeno.size();k++)
        {

            double allFamilyLikelihood=1.0;
            string fatherGenotype;
            string motherGenotype;


            for(int i=0;i<(int)current->spouse.size();i++)
            {

                Node* thisSpouse=current->spouse[i];
                double spouseLik=0;

                vector<string> SpouseGeno;
                if(thisSpouse->geno.compare("")==0)
                    SpouseGeno=geno;
                else
                    SpouseGeno.push_back(thisSpouse->geno);



                for(int l=0;l<(int)SpouseGeno.size();l++)
                {

                    double childLik=1.0;
                    for(int j=0;j<(int)current->child[thisSpouse].size();j++)
                    {

                        Node* thisChild=current->child[thisSpouse][j];


                        vector<string> ChildGeno;
                        if(thisChild->geno.compare("")==0)
                            ChildGeno=geno;
                        else
                            ChildGeno.push_back(thisChild->geno);



                        if(thisChild->ID.compare(from->ID)!=0)
                        {
                            double childGenoSum=0.0;
                            //cout<<"WELLLL\n";

                            for(int m=0;m<(int)ChildGeno.size();m++)
                            {
                                //childGenoSum+=(PhenoGivenGeno(thisChild->pheno,geno[m])*ChildGivenParent(geno[m],geno[l],geno[k])*RecursionLikelihood(thisChild,ped.first_founder,geno[m],"child"));

                                childGenoSum+=(PhenoGivenGeno(thisChild->pheno,ChildGeno[m],thisChild->sex,thisChild->father->geno,thisChild->mother->geno)*ChildGivenParent(ChildGeno[m],SpouseGeno[l],currentGeno[k])*RecursionLikelihood(thisChild,current,ChildGeno[m],"child"));

                            }
                            childLik*=childGenoSum;
                        }
                        else
                            childLik*=(ChildGivenParent(genotype,SpouseGeno[l],currentGeno[k]));


                    }
//                    if(current->ID.compare("34")==0)
//                        cout<<current->ID<<"\t"<<from->ID<<"\t"<<direction<<"\t"<<childLik<<"\t"<<genotype<<"\t"<<geno[l]<<"\t"<<geno[k]<<endl;

                fatherGenotype="";
                motherGenotype="";

                if(!thisSpouse->founder)
                {
                    fatherGenotype=thisSpouse->father->geno;
                    motherGenotype=thisSpouse->mother->geno;
                }
                spouseLik+=(childLik*PhenoGivenGeno(thisSpouse->pheno,SpouseGeno[l],thisSpouse->sex,fatherGenotype,motherGenotype)*RecursionLikelihood(thisSpouse,current,SpouseGeno[l],"spouse"));


                }

                allFamilyLikelihood*=spouseLik;

            }

            //finalLikelihood*=allSpouseLik;
            if(current->founder==false)
            {
                fatherGenotype=current->father->geno;
                motherGenotype=current->mother->geno;
                allFamilyLikelihood*=(PhenoGivenGeno(current->pheno,currentGeno[k],current->sex,fatherGenotype,motherGenotype)*RecursionLikelihood(current->father,current,currentGeno[k],"parent"));
            }
            else
            {


                allFamilyLikelihood*=(PhenoGivenGeno(current->pheno,currentGeno[k],current->sex,"","")*ProbGeno(currentGeno[k]));
            }

            finalLikelihood+=allFamilyLikelihood;
        }

         LikelihoodValue[current][genotype]=finalLikelihood;


        return finalLikelihood;
    }




    if(direction.compare("spouse")==0)
    {
        double allFamilyLik=1.0;
        string fatherGenotype="";
        string motherGenotype="";

        for(int i=0;i<(int)current->spouse.size();i++)
        {

            if(current->spouse[i]->ID.compare(from->ID)!=0)
            {
                Node* thisSpouse=current->spouse[i];
                double spouseLik=0;

                vector<string> SpouseGeno;
                if(thisSpouse->geno.compare("")==0)
                    SpouseGeno=geno;
                else
                    SpouseGeno.push_back(thisSpouse->geno);



                for(int l=0;l<(int)SpouseGeno.size();l++)
                {

                    double childLik=1.0;
                    for(int j=0;j<(int)current->child[thisSpouse].size();j++)
                    {

                        Node* thisChild=current->child[thisSpouse][j];

                        double childGenoSum=0.0;

                        vector<string> ChildGeno;
                        if(thisChild->geno.compare("")==0)
                            ChildGeno=geno;
                        else
                            ChildGeno.push_back(thisChild->geno);



                        for(int m=0;m<(int)ChildGeno.size();m++)
                        {
                            childGenoSum+=(PhenoGivenGeno(thisChild->pheno,ChildGeno[m],thisChild->sex,thisChild->father->geno,thisChild->mother->geno)*ChildGivenParent(ChildGeno[m],SpouseGeno[l],genotype)*RecursionLikelihood(thisChild,current,ChildGeno[m],"child"));
                        }
                        childLik*=childGenoSum;
                    }
                    fatherGenotype="";
                    motherGenotype="";

                    if(!thisSpouse->founder)
                    {
                        fatherGenotype=thisSpouse->father->geno;
                        motherGenotype=thisSpouse->mother->geno;
                    }
                    spouseLik+=(childLik*PhenoGivenGeno(thisSpouse->pheno,SpouseGeno[l],thisSpouse->sex,fatherGenotype,motherGenotype)*RecursionLikelihood(thisSpouse,current,SpouseGeno[l],"spouse"));
                }

            allFamilyLik*=spouseLik;
            }
        }


        if(current->founder==false)
        {
            allFamilyLik*=RecursionLikelihood(current->father,current,genotype,"parent");
           // cout<<current->ID<<"\t"<<genotype<<"\t"<<allFamilyLik<<endl;


        }


        else
            {
                //cout<<current->ID<<"\t"<<genotype<<"\t"<<allFamilyLik<<endl;
            allFamilyLik*=ProbGeno(genotype);
            }


            LikelihoodValue[current][genotype]=allFamilyLik;

return allFamilyLik;
    }




}


template <class Lik>
Lik Likelihood<Lik>::calLogLikelihood(Pedigree &ped)
{




    if(ped.first_founder==NULL)
        return 0.0;



    else
    {
        vector<string> founderGeno;
        if(ped.first_founder->geno.compare("")==0)
            founderGeno=geno;
        else
            founderGeno.push_back(ped.first_founder->geno);


        double finalLikelihood=0;
        for(int k=0;k<(int)founderGeno.size();k++)
        {
            double allSpouseLik=1;
            for(int i=0;i<(int)ped.first_founder->spouse.size();i++)
            {

                Node* thisSpouse=ped.first_founder->spouse[i];
                double spouseLik=0;

                vector<string> SpouseGeno;
                if(thisSpouse->geno.compare("")==0)
                    SpouseGeno=geno;
                else
                    SpouseGeno.push_back(thisSpouse->geno);

                for(int l=0;l<(int)SpouseGeno.size();l++)
                {

                    double childLik=1;
                    //RecursionLikelihood()
                    //countSamples(first_founder->spouse[i],first_founder->ID,"spouse",count);
                    for(int j=0;j<(int)ped.first_founder->child[thisSpouse].size();j++)
                    {

                        Node* thisChild=ped.first_founder->child[thisSpouse][j];
                        double childGenoSum=0;

                        vector<string> ChildGeno;
                        if(thisChild->geno.compare("")==0)
                            ChildGeno=geno;
                        else
                            ChildGeno.push_back(thisChild->geno);


                        for(int m=0;m<(int)ChildGeno.size();m++)
                        {
                            //cout<<"WELLLLLLLL\n";

                            childGenoSum+=(PhenoGivenGeno(thisChild->pheno,ChildGeno[m],thisChild->sex,thisChild->father->geno,thisChild->mother->geno)*ChildGivenParent(ChildGeno[m],SpouseGeno[l],founderGeno[k])*RecursionLikelihood(thisChild,ped.first_founder,ChildGeno[m],"child"));
                            //cout<<"wat\n";
                        }
                         childLik*=childGenoSum;
                    }
                    string fatherGenotype="";
                    string motherGenotype="";

                    if(!thisSpouse->founder)
                    {
                        fatherGenotype=thisSpouse->father->geno;
                        motherGenotype=thisSpouse->mother->geno;
                    }
                    spouseLik+=(childLik*PhenoGivenGeno(thisSpouse->pheno,SpouseGeno[l],thisSpouse->sex,fatherGenotype,motherGenotype)*RecursionLikelihood(thisSpouse,ped.first_founder,SpouseGeno[l],"spouse"));


                }
                allSpouseLik*=spouseLik;
            }

            //cout<<ped.first_founder->ID<<"\t"<<geno[k]<<"\t"<<allSpouseLik<<endl;
            finalLikelihood+=(PhenoGivenGeno(ped.first_founder->pheno,founderGeno[k],ped.first_founder->sex,"","")*allSpouseLik*ProbGeno(founderGeno[k]));
        }


    if(finalLikelihood<ZEPS)
        return log(ZEPS);
    else
        return log(finalLikelihood);

    }

}




#endif
