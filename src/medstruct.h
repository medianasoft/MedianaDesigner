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

struct MeanSD 
{
    double mean;
    double sd;    
};

struct ModelInformation {
    int model_index;
    int n_parameters;
    vector<double> initial_values;
    vector<double> coef;
    int status;
    double criterion;
    double target_dose;
    double convergence_criterion;
    vector<double> bounds;
};

struct AllSurvivalData {
    vector<int> stratum;
    vector<double> start;
    vector<double> dropout;
    vector<double> os;
    vector<double> os_local;
    vector<double> os_local_censor;
    vector<double> os_local_start;
};

#endif // MEDDESSTRUCT_H
