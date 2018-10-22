#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
    return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_MATRIX(Pounds);
    DATA_MATRIX(IOA);
    DATA_INTEGER(n);
    DATA_INTEGER(NRegion);

    PARAMETER_VECTOR(logq);
    PARAMETER(logk);
    PARAMETER(logr);
    PARAMETER(logSigmaO);
    PARAMETER(dummy);

    vector<Type> q = exp(logq);
    Type k = exp(logk);
    Type r = exp(logr);
    Type SigmaO = exp(logSigmaO);
    
    matrix<Type> B(n+1,NRegion);
    matrix<Type> Ifit(n,NRegion);
    //vector<Type> B(n+1);
    //vector<Type> Ifit(n);
    Type obj_fun = 0.0;
    
    B(0,0) = k;
    B(0,1) = k;
    
    for (int region = 0; region < NRegion; region++) {
        for (int t = 0; t < n; t++) {   // r is region, t is year
            B(t+1, region) = B(t,region) + r*B(t,region)*(1-(B(t,region)/k))-Pounds(t,region);
            Ifit(t,region) = q(region)*B(t,region);
        }
     
    }

    vector<Type> temp(n);
    for(int region =0; region <NRegion; region++){
    temp = 0;
    
    for (int t=0; t<n; t++) {
        if(!isNA(IOA(t))){
        
            temp(region) += -dnorm(log(Ifit(t)), log(IOA(t)), SigmaO, true);
            
        }
        }
        obj_fun += temp(region);
    }
   
    REPORT(B);
    REPORT(Ifit);
    REPORT(logq);
    REPORT(SigmaO);
    REPORT(q);
    REPORT(obj_fun);
    return obj_fun;

}
