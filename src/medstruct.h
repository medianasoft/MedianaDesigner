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

struct AllSurvivalData {
    vector<int> stratum;
    vector<double> start;
    vector<double> dropout;
    vector<double> os;
    vector<double> os_local;
    vector<double> os_local_censor;
    vector<double> os_local_start;
};

struct MeanSD 
{
    double mean;
    double sd;    
};

// Information required to fit a dose-response model and fitting results
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


/*
//From ADRand


// Crossing probabilities
struct CrossProb
{
    vector<double> problo;
    vector<double> probhi;    
};

struct GSDataModel
{
    int endpoint_distribution;
    int endpoint_test;
    vector<int> sample_size;
    vector<int> event_count;
    vector<double> control_parameter; 
    vector<double> treatment_parameter;
    double enrollment_period;
    int enrollment_distribution;
    vector<double> enrollment_parameter;
    int dropout_distribution;
    vector<double> dropout_parameter;
    double followup_period;
    int nsims;
    int time_to_event;
};

struct GSAnalysisModel
{
    vector<double> control_parameter; 
    vector<double> treatment_parameter;
    int endpoint_distribution;
    int endpoint_test;
    int direction;
    vector<double> info_frac;
    int upper_spending_family;    
    double upper_spending_parameter;    
    int lower_spending_family;    
    double lower_spending_parameter;    
    double alpha;
    double beta;
};

struct GSResults
{
    vector<double> upper_stat;
    vector<double> upper_p;
    vector<double> upper_es;
    vector<double> upper_hr;
    vector<double> lower_stat;
    vector<double> lower_p;
    vector<double> lower_es;
    vector<double> lower_hr;
    vector<double> upper_stop_prob_h1;
    vector<double> upper_stop_prob_h0;
    vector<double> lower_stop_prob_h1;
    vector<double> lower_stop_prob_h0;
    vector<double> gs_sample;
    vector<double> sample_size;
    vector<double> sim_look_time;
    vector<double> sim_upper_stop_prob;
    vector<double> sim_lower_stop_prob;
    vector<double> sim_nenrolled;
    vector<double> sim_ncompleters;
    double sim_nexpenrolled;
    double sim_nexpcompleters;
    vector<double> sim_noverrun;
    vector<double> sim_ndropout;

};



// Parameters of an adaptive design with population selection
struct ADPSParameters
{
    int endpoint_distribution;
    int direction;    
    int endpoint_test;
    int time_to_event;
    vector<double> control_parameters_plus;
    vector<double> treatment_parameters_plus;
    vector<double> control_parameters_minus;
    vector<double> treatment_parameters_minus;
    double alpha;
    double beta;
    int alpha_spending_family;
    double alpha_spending_parameter;
    int beta_spending_family;
    double beta_spending_parameter;
    double enrollment_period;
    double followup_period;
    int enrollment_distribution;
    int dropout_distribution;
    vector<double> enrollment_parameter;
    vector<double> dropout_parameter;
    vector<int> sample_size;
    vector<int> event_count;
    vector<int> n;
    int nsims;
    int objective;
    double influence_threshold;
    double interaction_threshold;
    vector<int> design_flag;
    vector<double> DesignA;
    vector<double> DesignB;
    double precision;

};

// Simulation results for an adaptive design with population selection
struct ADPSResults
{
    vector<double> common_design_characteristics;
    vector<double> common_design_parameters;
    vector<double> DesignA_characteristics;
    vector<double> DesignB_characteristics;

};


*/
#endif // MEDDESSTRUCT_H
