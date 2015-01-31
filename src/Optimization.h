



#ifndef __OPTIMIZATION__
#define __OPTIMIZATION__
#include <vector>
#include <cmath>
#include <iostream>
#define EPS 1e-10

// Simplex contains (dim+1)*dim points

template <class F>
class Optimization
{
    protected:
        std::vector<std::vector<double> > X;
        std::vector<double> Y;
        std::vector<double> midPoint;
        std::vector<double> thruLine;
        int iter; // No. of Iterations
        int dim, idxLo, idxHi, idxNextHi;
        void evaluateFunction(F& LikeFunction);
        void evaluateExtremes();
        void prepareUpdate();
        bool updateSimplex(F& LikeFunction, double scale);
        void contractSimplex(F& LikeFunction);
        static int check_tol(double fmax, double fmin, double ftol);
    public:
        int iterations()
        {
            return iter;
        }
        Optimization(double* p, int d);
        void amoeba(F& LikeFunction, double tol);
        std::vector<double>& xmin();
        double ymin();
};


template <class F>
Optimization<F>::Optimization(double* p, int d) : dim(d)
{
    iter=0;
    X.resize(dim+1);
    Y.resize(dim+1);
    midPoint.resize(dim);
    thruLine.resize(dim);
    for(int i=0; i < dim+1; ++i)
        X[i].resize(dim);

    for(int i=0; i < dim+1; ++i)
        for(int j=0; j < dim; ++j)
            X[i][j] = p[j];

    for(int i=0; i < dim; ++i)
        X[i][i] += 1.;
}

template <class F>
void Optimization<F>::evaluateFunction(F& LikeFunction)
{
    for(int i=0; i < dim+1; ++i)
        Y[i] = LikeFunction(X[i]);

}


template <class F>
void Optimization<F>::evaluateExtremes()
{
    if ( Y[0] > Y[1] )
    {
        idxHi = 0;
        idxLo = idxNextHi = 1;
    }
    else
    {
        idxHi = 1;
        idxLo = idxNextHi = 0;
    }

    for(int i=2; i < dim+1; ++i)
    {
        if ( Y[i] <= Y[idxLo] )
            idxLo = i;
        else if ( Y[i] > Y[idxHi] )
        {
            idxNextHi = idxHi;
            idxHi = i;

        }
        else if ( Y[i] > Y[idxNextHi] )
            idxNextHi = i;

    }
}


template <class F>
void Optimization<F>::prepareUpdate()
{
    for(int j=0; j < dim; ++j)
        midPoint[j] = 0;

    for(int i=0; i < dim+1; ++i)
        if ( i != idxHi )
            for(int j=0; j < dim; ++j)
                midPoint[j] += X[i][j];

    for(int j=0; j < dim; ++j)
    {
        midPoint[j] /= dim;
        thruLine[j] = X[idxHi][j] - midPoint[j];
    }
}


template <class F>
bool Optimization<F>::updateSimplex(F& LikeFunction, double scale)
{
    std::vector<double> nextPoint;
    nextPoint.resize(dim);
    for(int i=0; i < dim; ++i)
        nextPoint[i] = midPoint[i] + scale * thruLine[i];

    double fNext = LikeFunction(nextPoint);
    if ( fNext < Y[idxHi] )
    { // exchange with maximum
        for(int i=0; i < dim; ++i)
            X[idxHi][i] = nextPoint[i];

        Y[idxHi] = fNext;
        return true;
    }
    else
        return false;
}


template <class F>
void Optimization<F>::contractSimplex(F& LikeFunction)
{
    for(int i=0; i < dim+1; ++i)
    {
    if ( i != idxLo )
    {
        for(int j=0; j < dim; ++j)
            X[i][j] = 0.5*( X[idxLo][j] + X[i][j] );
        Y[i] = LikeFunction(X[i]);
        }
    }
}




template <class F>
void Optimization<F>::amoeba(F& LikeFunction, double tol)
{
    evaluateFunction(LikeFunction);


    while(true)
    {
        iter++;
        evaluateExtremes();
        prepareUpdate();
        if (check_tol(Y[idxHi],Y[idxLo],tol) ) break;
        updateSimplex(LikeFunction, -1.0); // reflection
        //cout<<Y[idxLo]<<endl;
        if ( Y[idxHi] < Y[idxLo] )
        {
            updateSimplex(LikeFunction, -2.0); // expansion
        }
        else if ( Y[idxHi] >= Y[idxNextHi] )
            if ( !updateSimplex(LikeFunction, 0.5) )
                contractSimplex(LikeFunction);
    }
}



template <class F>
std::vector<double>& Optimization<F>::xmin()
{
    return X[idxLo];
}


template <class F>
double Optimization<F>::ymin()
{
    return Y[idxLo];
}


template <class F>
int Optimization<F>::check_tol(double fmax, double fmin, double ftol)
{
    double delta = fabs(fmax - fmin); double accuracy = (fabs(fmax) + fabs(fmin)) * ftol;
    return (delta < (accuracy + EPS));
}
#endif // __OPTIMIZATION__

