#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_VECTOR(Pounds);
    DATA_VECTOR(IOA);
    int n = Pounds.size();

    PARAMETER(logq);
    PARAMETER(logk);
    PARAMETER(logr);
    PARAMETER(logSigmaO);

    Type q = exp(logq);
    Type k = exp(logk);
    Type r = exp(logr);
    Type SigmaO = exp(logSigmaO);
    
    //int n1 = 0;   //local integer
    //n1 = n + 1;
    vector<Type> B(n+1);
    vector<Type> Ifit(n);
    Type obj_fun;
    
    B(0) = k;
    for (int t = 0; t < n; t++) {
        B(t+1) = B(t) + r*B(t)*(1-(B(t)/k))-Pounds(t);
        Ifit(t) = q*B(t);
    }

    
    obj_fun = -sum(dnorm(log(IOA), log(Ifit), SigmaO, true));

    REPORT(B);
    REPORT(Ifit);
    
    return obj_fun;

}
