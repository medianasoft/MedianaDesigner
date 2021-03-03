// All data structures used in the package 

#ifndef MEDSTRUCT_H
#define MEDSTRUCT_H

struct OutcomeCensor
{
    vector<double> outcome;
    vector<double> censor;    
};

struct EfficacyFutility
{
    vector<double> efficacy;
    vector<double> futility;    
};

struct TestResult 
{
    double pvalue;
    double test_stat;    
};


#endif // MEDDESSTRUCT_H
