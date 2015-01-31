
#ifndef __NUMERICAL_H
#define __NUMERICAL_H

#include <ctime>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#define REPS 1e-30
using namespace std;



template <class H>
class PartialDifferentiation
{

    private:

        double diff;
        int countCheck;
    public:

        int type; // 1: d2/dx2, 2: d/dx(d/dy), 3: d/dy(d/dx), 4L d2/dy2
        H func;
          double relError;
        PartialDifferentiation()
        {

        }
          double varPar;
        PartialDifferentiation(H &function,   double relE,  double var,int t)
        {
            countCheck=0;
            func=function;
            relError=relE/10;
            varPar=var;
            type=t;

        }
          double operator [](  double x);



        bool checktol(  double difference);
};



template <class F>
class Differentiation
{

    private:

          double diff;
    public:
        F func;
          double relError;
        Differentiation()
        {

        }
        Differentiation(F &functions,   double relE)
        {
            func=functions;
            relError=relE;
        };
          double operator [](  double x);
        double SecondDiff(double x)
        {

            diff=1e-2;
            double f0=func[x];
            double fminh=func[x-diff];
            double fpluh=func[x+diff];

            double deriv=0.0;
            double deriv_new=0.0;
            deriv=(fpluh-(2*f0)+(fminh))/(diff*diff);

            while(diff>REPS)
            {
                diff/=2;
                fminh=func[x-diff];
                fpluh=func[x+diff];

                deriv_new=(fpluh-(2*f0)+(fminh))/(diff*diff);
        //cout<<x+diff<<"\tf\'("<<x+diff<<") = "<<temp<<"\t"<<x<<"\tf\'("<<x<<") = "<<temp2<<"\t"<<deriv_new<<"\t"<<deriv<<"\t"<<deriv_new-deriv<<endl;

        //cout<<deriv<<" "<<deriv_new<<endl;
        if(checktol(deriv_new-deriv)==true)
            {
                //cout<<endl;
                return deriv_new;
            }
        deriv=deriv_new;

    }
    return deriv_new;


        };



        double SecondDiffMult(vector<double> xx,int type)
        {

            vector<  double> par0(2),parph(2),parmh(2);
            diff=1e-2;
            int countCheck=0;
            double deriv=0.0;
            double deriv_new=0.0;
            par0=xx;
            parph=xx;
            parmh=xx;

            if(type==1)
            {
                parmh[0]-=diff;
                parph[0]+=diff;

            }
            else if(type==4)
            {
                        parmh[1]-=diff;
                        parph[1]+=diff;
            }
            deriv=(func(parph)-(2*func(par0))+(func(parmh)))/(diff*diff);
            while(diff>REPS)

            {
                diff/=2;
                if(type==1)
                {
                    parph[0]=xx[0]+diff;
                    parmh[0]=xx[0]-diff;
                }
                else if(type==4)
               {
                    parph[1]=xx[1]+diff;
                    parmh[1]=xx[1]-diff;
                }

                deriv_new=(func(parph)-(2*func(par0))+(func(parmh)))/(diff*diff);
                if(checktol(deriv_new-deriv)==true)
                    {
                        countCheck++;
                        if(countCheck==8)
                        return deriv_new;
                    }
                else
                    countCheck=0;

            deriv=deriv_new;

    }

    return deriv_new;
}





        bool checktol(  double difference);
};


template <class G>
class SecondDifferentiation
{

    private:

    public:
        SecondDifferentiation()
        {

        }
        vector<vector<  double> > doubleDifferential(G &func,vector<  double> x,  double relError);
};


template <class G>
vector<vector<  double> >  SecondDifferentiation<G>::doubleDifferential(G &func,vector<  double> x,  double relError)
{


    vector<vector<   double> > SecondDeriv(x.size());


    if(x.size()==1)
    {

        Differentiation<G> FirstDiff(func,relError);
        SecondDeriv[0].push_back(FirstDiff.SecondDiff(x[0]));
//        Differentiation<Differentiation<G> > SecondDiff(FirstDiff,relError);
//        SecondDeriv[0].push_back(SecondDiff[x[0]]);
        return SecondDeriv;
    }


    for(int i=0;i<(int)x.size();i++)
    {
        SecondDeriv[i].resize(x.size());
        for(int j=0;j<(int)x.size();j++)
        {

            if(i==j)
            {
                if(i==0)
                {
                     Differentiation<G> FirstDiff(func,relError);
                     SecondDeriv[i][j]=FirstDiff.SecondDiffMult(x,1);
                }
                else
                {
                    Differentiation<G> FirstDiff(func,relError);
                     SecondDeriv[i][j]=FirstDiff.SecondDiffMult(x,4);
                }

            }
            else
            {

                if(i==0)
                {
                    PartialDifferentiation<G> Partial(func,relError,x[1],2);
                    Differentiation<PartialDifferentiation<G> > SecondPartialDiff(Partial,relError);
                    SecondDeriv[i][j]=SecondPartialDiff[x[0]];
                }
                else
                {
                    PartialDifferentiation<G> Partial(func,relError,x[0],3);
                    Differentiation<PartialDifferentiation<G> > SecondPartialDiff(Partial,relError);
                    SecondDeriv[i][j]=SecondPartialDiff[x[1]];
                    }
            }
        }
    }

//    SecondDeriv[0][1]=(SecondDeriv[0][1]+SecondDeriv[1][0])/2;
//    SecondDeriv[1][0]=SecondDeriv[0][1];

    for(int i=0;i<(int)x.size();i++)
    {
        for(int j=0;j<(int)x.size();j++)
            SecondDeriv[i][j]=floor(SecondDeriv[i][j]*10000+0.5)/10000;
    }

    return SecondDeriv;



}





template <class H>
bool PartialDifferentiation<H>::checktol(  double difference)
{

    if(abs(difference)<relError)
        return true;
    else
        return false;

}


template <class H>
  double PartialDifferentiation<H>::operator [](double x)
{
    int countCheck=0;
    vector<  double> par0(2),par1(2),parph(2),parmh(2),par2ph(2),par2mh(2);
    diff=1e-2;
    double deriv=0.0;
    double deriv_new=0.0;


    if(type==2)
    {
        par0[0]=x;
        par0[1]=varPar;

        parph[0]=x;
        parph[1]=varPar+diff;
        parmh[0]=x;
        parmh[1]=varPar-diff;

        par2ph[0]=x;
        par2ph[1]=varPar+diff+diff;
        par2mh[0]=x;
        par2mh[1]=varPar-diff-diff;
    }
    else if(type==3)
    {
        par0[1]=x;
        par0[0]=varPar;

        parph[1]=x;
        parph[0]=varPar+diff;
        parmh[1]=x;
        parmh[0]=varPar-diff;

        par2ph[1]=x;
        par2ph[0]=varPar+diff+diff;
        par2mh[1]=x;
        par2mh[0]=varPar-diff-diff;
    }



    deriv=((func(par2mh))-(8*func(parmh))+(8*func(parph))-(func(par2ph)))/(12*diff);

    while(diff>REPS)
    {
        diff/=2;
        if(type==2)
        {
            parph[1]=varPar+diff;
            parmh[1]=varPar-diff;
            par2ph[1]=varPar+diff+diff;
            par2mh[1]=varPar-diff-diff;
        }
        else if(type==3)
        {
            parph[0]=varPar+diff;
            parmh[0]=varPar-diff;
            par2ph[0]=varPar+diff+diff;
            par2mh[0]=varPar-diff-diff;
        }

        deriv_new=((func(par2mh))-(8*func(parmh))+(8*func(parph))-(func(par2ph)))/(12*diff);

        if(checktol(deriv_new-deriv)==true)
        {
            countCheck++;
            if(countCheck==3)
                return deriv_new;
        }
        else
            countCheck=0;

        deriv=deriv_new;

    }

    return deriv_new;
}



template <class F>
  double Differentiation<F>::operator [](  double x)
{

    int countCheck=0;
    diff=1e-2;
    double deriv=0.0;
    double deriv_new=0.0;
    deriv=(func[x+diff]-func[x-diff])/(2*diff);
    //cout<<deriv<<endl;
    while(diff>REPS)
    {
        diff/=2;
        double temp=func[x+diff],temp2=func[x-diff];
        deriv_new=(temp-temp2)/(2*diff);
        //cout<<deriv_new<<endl;
        if(checktol(deriv_new-deriv)==true)
        {
            countCheck++;
            if(countCheck==3)
                return deriv_new;
        }
        else
            countCheck=0;
        deriv=deriv_new;

    }
    return deriv_new;
}




template <class F>
bool Differentiation<F>::checktol(  double difference)
{

    if(abs(difference)<relError)
        return true;
    else
        return false;

}
#endif


