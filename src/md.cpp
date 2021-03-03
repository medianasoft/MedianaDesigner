#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

#include "medstruct.h"
#include "medsupport.h"
#include "meddist.h"
#include "medstattest.h"
#include "meddebug.h"

#include <exception>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <math.h> 
#include <iomanip>
#include <fstream>
#include <string>

struct AllSurvivalData {
    vector<int> stratum;
    vector<double> start;
    vector<double> dropout;
    vector<double> os;
    vector<double> os_local;
    vector<double> os_local_censor;
    vector<double> os_local_start;
};

vector<double> HochbergOutcome(const vector<double> &pvalue, const double &alpha) {

    int m = pvalue.size();
    vector<double> outcome(2);

    if (m == 2) {

        if (pvalue[0] <= pvalue[1]) {

            if (pvalue[1] <= alpha) {
                outcome[0] = 1;
                outcome[1] = 1;
            }

            if (pvalue[1] > alpha && pvalue[0] <= alpha / 2) { 
                outcome[0] = 1;
                outcome[1] = 0;
            }

        }

        if (pvalue[0] > pvalue[1]) {

            if (pvalue[0] <= alpha) {
                outcome[0] = 1;
                outcome[1] = 1;
            }

            if (pvalue[0] > alpha && pvalue[1] <= alpha / 2) { 
                outcome[0] = 0;
                outcome[1] = 1;
            }

        }

    }

    return outcome;

}

struct MeanSD 
{
    double mean;
    double sd;    
};

double ComputeRate(const vector<double> &vec) {

    double sum = std::accumulate(vec.begin(), vec.end(), 0);
    return sum / (vec.size() + 0.0);
}

MeanSD ComputeMeanSD(const std::vector<double> &x) {

    double Sm=0.0;       // Sm1 = Sample 1 Mean.
    double Var=0.0;       // Var1 = Sample 1  Variation.
    unsigned Sn = x.size();      // Sn1 = Sample 1 Size.
    int i;
    double j;
    for (i = 0; i<Sn ; ++i)
    {
        j = x[i];
        Sm += j;
        Var += j*j;
    }

    MeanSD res;

    Sm=Sm/Sn;
    res.mean = Sm;
    res.sd = sqrt((Var/Sn-Sm*Sm)*Sn/(Sn-1));

    return res;

}

double ComputeEffectSize(const std::vector<double> &sample1, const std::vector<double> &sample2, const int &endpoint_distribution, const int &endpoint_test) {

    double effect_size = 0.0;
    double n1 = sample1.size();
    double n2 = sample2.size();
    double pooled_sd;

    // Normal distribution
    if (endpoint_distribution == 1) {

        MeanSD MeanSD1 = ComputeMeanSD(sample1);
        MeanSD MeanSD2 = ComputeMeanSD(sample2);

        pooled_sd = sqrt(((n1-1.0) * MeanSD1.sd * MeanSD1.sd  + (n2-1.0) * MeanSD2.sd * MeanSD2.sd) / (n1 + n2 - 2.0));

        effect_size = (MeanSD2.mean - MeanSD1.mean) / pooled_sd;

    }

    // Binary distribution
    if (endpoint_distribution == 2) {

        double estimatex = ComputeRate(sample1);
        double estimatey = ComputeRate(sample2);
        double ave = (estimatex+estimatey)/2.0;

        // Without pooled variance
        if (endpoint_test == 1) effect_size = abs(estimatey - estimatex) / sqrt(estimatex*(1.0-estimatex)+ estimatey*(1.0-estimatey)); 

        // With pooled variance
        if (endpoint_test == 2) effect_size = abs(estimatey - estimatex) / sqrt(ave*(1.0-ave)); 

    }


    return effect_size;

} 

// Extract individual samples from overall sample
vector<double> ExtractSamples(const vector<double> overall_sample, const int &start, const int &end) {

    vector<double> result;
    
    result.clear();
    result.insert(result.end(), overall_sample.begin() + start, overall_sample.begin() + end);

    return result;

}

vector<double> CombineVec(const vector<double> &x, const vector<double> &y) {

    vector<double> z = x;

    z.insert(z.end(), y.begin(), y.end());

    return z;

}

// Direction 1: higher values are beneficial, 2: lower values are beneficial 
double CondPower(const double &test_stat, const int &n_interim, const int &n_final, const int &direction, const int &observed, const double &assumed_effect_size, const double &alpha) {

    double z_alpha = rcpp_qnorm(1.0 - alpha);

    // Derive the effect size from test_stat
    double effect_size;

    if (observed == 1) effect_size = test_stat / sqrt(n_interim / 4.0);
    if (observed == 0) effect_size = assumed_effect_size;

    double w1 = (n_interim + 0.0) / (n_final + 0.0);
    double w2 = 1.0 - w1;

    double val = effect_size * sqrt((n_final - n_interim) / 4.0) + test_stat * sqrt(w1 / w2) - z_alpha / sqrt(w2);

    return rcpp_pnorm(val); 

}

// Compute the total updated event count (increase after the interim analysis) for the standard adaptive design
int UpdatedEventCount(const double &test_stat, const int &n_interim, const int &n_final, const int &n_maximum, const int &direction, const int &observed, const double &assumed_effect_size, const double &target_power, const double &alpha) {

    double z_alpha = rcpp_qnorm(1.0 - alpha);
    double z_beta = rcpp_qnorm(target_power);

    double effect_size;

    if (observed == 1) effect_size = test_stat / sqrt(n_interim / 4.0);
    if (observed == 0) effect_size = assumed_effect_size;
   
    double w1 = (n_interim + 0.0) / (n_final + 0.0);
    double w2 = 1.0 - w1;

    double temp, n_updated;

    int n;

    if (abs(effect_size) > 0.00001) {

        temp = (z_alpha / sqrt(w2) - test_stat * sqrt(w1 / w2) + z_beta) / effect_size;
        n_updated = 4.0 * temp * temp;
        n = (int) (n_updated + 1); 

    } 
    else {

        n = 0;

    }

    if (n > n_maximum - n_interim) n = n_maximum - n_interim;
    if (n < n_final - n_interim) n = n_final - n_interim;

    return n;

}

// Apply hypothesis selection rules based on the influence and interaction conditions in adaptive population selection trials
vector<double> HypothesisSelection(const double &effect_size_minus, const double &effect_size_plus, const double &influence_threshold, const double &interaction_threshold) {

    int influence_flag = 0, interaction_flag = 0;

    double outcome1 = 0.0, outcome2 = 0.0, outcome3 = 0.0;

    vector<double> res(3);

    if (effect_size_minus >= influence_threshold) influence_flag = 1;

    if (interaction_threshold >= 0) {

        // Interaction condition
        if (influence_flag == 1) {
            if (effect_size_plus / effect_size_minus >= interaction_threshold) interaction_flag = 1;
        }                

        if (influence_flag == 1 && interaction_flag == 0) outcome1 = 1.0;

        if (influence_flag == 0) outcome2 = 1.0;

        if (influence_flag == 1 && interaction_flag == 1) outcome3 = 1.0;

    } else {

        // Both populations are always selected if the interaction condition is not specified
        outcome3 = 1.0;

    }

    // Save the results

    // Final analysis will be performed in the overall population only 
    res[0] = outcome1;

    // Final analysis will be performed in the biomarker-positive subpopulation only 
    res[1] = outcome2;

    // Final analysis will be performed in both populations 
    res[2] = outcome3;

    return res;

}
// End of HypothesisSelection

// [[Rcpp::export]]
List ADSSModC(const List &parameters_arg) {

    List parameters(parameters_arg);

    vector<int> sample_size = as<vector<int>>(parameters["sample_size"]);
    vector<int> sample_size_ia1 = as<vector<int>>(parameters["sample_size_ia1"]);
    vector<int> sample_size_ia2 = as<vector<int>>(parameters["sample_size_ia2"]);
    vector<int> sample_size_fa1 = as<vector<int>>(parameters["sample_size_fa1"]);
    vector<int> sample_size_fa2 = as<vector<int>>(parameters["sample_size_fa2"]);

    int event_count_ia1 = as<int>(parameters["event_count_ia1"]);
    int event_count_ia2 = as<int>(parameters["event_count_ia2"]);
    int event_count_fa1 = as<int>(parameters["event_count_fa1"]);
    int event_count_fa2 = as<int>(parameters["event_count_fa2"]);

    vector<double> dropout_parameter = as<vector<double>>(parameters["dropout_parameter"]);

    int endpoint_index = as<int>(parameters["endpoint_index"]);

    vector<double> means = as<vector<double>>(parameters["means"]);
    vector<double> sds = as<vector<double>>(parameters["sds"]);
    vector<double> rates = as<vector<double>>(parameters["rates"]);
    vector<double> hazard_rates = as<vector<double>>(parameters["hazard_rates"]);

    double futility_threshold = as<double>(parameters["futility_threshold"]);
    vector<double> underpowered_zone = as<vector<double>>(parameters["promising_interval"]);
    double weight = as<double>(parameters["weight"]);

    double target_power = as<double>(parameters["target_power"]);

    double enrollment_period = as<double>(parameters["enrollment_period"]);
    int enrollment_distribution = as<int>(parameters["enrollment_distribution"]);
    double enrollment_parameter = as<double>(parameters["enrollment_parameter"]);

    int nsims = as<int>(parameters["nsims"]);
    double alpha = as<double>(parameters["alpha"]);

    /*******************************************************************/

    double stage1_test, stage2_test, trad_outcome, ratio = (double) sample_size[0] / (double(sample_size[0] + sample_size[1])), pvalue1, pvalue2;

    int i, sim, current_stratum, ntotal = sample_size[0] + sample_size[1], sample_size_increase;

    AllSurvivalData survival_data;

    OutcomeCensor outcome_censor, outcome_censor_control, outcome_censor_treatment;

    TestResult test_result1, test_result2;

    vector<double> look_time(3), outcome(4), event_count(3), dropout, hr(3), cp(2), overall_control_sample, overall_treatment_sample, control_sample, treatment_sample;

    vector<int> stratum_list, sample_list;

    NumericMatrix sim_results(nsims, 19);

    /*******************************************************************/

    for (sim = 0; sim < nsims; sim++) {

        hr = fillvec(3, 0.0);
        cp = fillvec(2, 0.0);

        look_time = fillvec(3, 0.0);
        outcome = fillvec(4, 0.0); 
        event_count = fillvec(3, 0.0);

        pvalue1 = 0.0;
        pvalue2 = 0.0;
        sample_size_increase = 0;

        // Normal and binary endpoints
        if (endpoint_index == 1 || endpoint_index == 2) {

            // Control and treatment samples
            if (endpoint_index == 1) {
                overall_control_sample = Normal(sample_size_fa2[0], means[0], sds[0]);
                overall_treatment_sample = Normal(sample_size_fa2[1], means[1], sds[1]);
            }
            if (endpoint_index == 2) {
                overall_control_sample = Binary(sample_size_fa2[0], rates[0]);
                overall_treatment_sample = Binary(sample_size_fa2[1], rates[1]);
            }

            // Outcomes
            // 0: Futility at IA1
            // 1: Significance at IA2
            // 2: Increase the number of events
            // 3: Significance at FA

            // Traditional analysis

            sample_list.clear();
            sample_list.push_back(1);        
            sample_list.push_back(2);    
            sample_list.push_back(3); 

            control_sample = ExtractSamples(overall_control_sample, 0, sample_size_fa1[0]);
            treatment_sample = ExtractSamples(overall_treatment_sample, 0, sample_size_fa1[1]);

            if (endpoint_index == 1) test_result1 = TTest(control_sample, treatment_sample, 0.0, 1);
            if (endpoint_index == 2) test_result1 = PropTest(control_sample, treatment_sample, 0.0, 1);

            pvalue1 = test_result1.pvalue;

            trad_outcome = (pvalue1 <= alpha);        

            // Interim analysis 1      

            control_sample = ExtractSamples(overall_control_sample, 0, sample_size_ia1[0]);
            treatment_sample = ExtractSamples(overall_treatment_sample, 0, sample_size_ia1[1]);

            if (endpoint_index == 1) test_result1 = TTest(control_sample, treatment_sample, 0.0, 1);
            if (endpoint_index == 2) test_result1 = PropTest(control_sample, treatment_sample, 0.0, 1);

            // Conditional power
            cp[0] = CondPower(test_result1.test_stat, sample_size_ia1[0] + sample_size_ia1[1], sample_size_fa1[0] + sample_size_fa1[1], 1, 1, 0.0, alpha);

            outcome[0] = (cp[0] <= futility_threshold);

            // Interim analysis 2

            if (outcome[0] < 1.0) {

                control_sample = ExtractSamples(overall_control_sample, 0, sample_size_ia2[0]);
                treatment_sample = ExtractSamples(overall_treatment_sample, 0, sample_size_ia2[1]);

                if (endpoint_index == 1) test_result1 = TTest(control_sample, treatment_sample, 0.0, 1);
                if (endpoint_index == 2) test_result1 = PropTest(control_sample, treatment_sample, 0.0, 1);

                // Conditional power
                cp[1] = CondPower(test_result1.test_stat, sample_size_ia2[0] + sample_size_ia2[1], sample_size_fa1[0] + sample_size_fa1[1], 1, 1, 0.0, alpha);

                pvalue1 = test_result1.pvalue;

            }

            if (outcome[0] < 1.0 && outcome[1] < 1.0) outcome[2] = (cp[1] > underpowered_zone[0] && cp[1] <= underpowered_zone[1]); 

            // Final analysis (before sample size adjustment)
            if (outcome[0] < 1.0 && outcome[1] < 1.0 && outcome[2] < 1.0) {

                control_sample = ExtractSamples(overall_control_sample, 0, sample_size_fa1[0]);
                treatment_sample = ExtractSamples(overall_treatment_sample, 0, sample_size_fa1[1]);

                if (endpoint_index == 1) test_result2 = TTest(control_sample, treatment_sample, 0.0, 1);
                if (endpoint_index == 2) test_result2 = PropTest(control_sample, treatment_sample, 0.0, 1);

                pvalue2 = test_result2.pvalue;

                outcome[3] = (pvalue2 <= alpha);        

            } 

            // Final analysis (after sample size adjustment)
            if (outcome[0] < 1.0 && outcome[1] < 1.0 && outcome[2] > 0.0) {

                // Updated sample size in both arms
                sample_size_increase = UpdatedEventCount(test_result1.test_stat, sample_size_ia2[0] + sample_size_ia2[1], sample_size_fa1[0] + sample_size_fa1[1], sample_size_fa2[0] + sample_size_fa2[1], 1, 1, 0.0, target_power, alpha); 

                control_sample = ExtractSamples(overall_control_sample, sample_size_ia2[0], sample_size_ia2[0] + (int) (ratio * sample_size_increase));
                treatment_sample = ExtractSamples(overall_treatment_sample, sample_size_ia2[1], sample_size_ia2[1] + (int) ((1.0 - ratio) * sample_size_increase));

                if (endpoint_index == 1) test_result2 = TTest(control_sample, treatment_sample, 0.0, 1);
                if (endpoint_index == 2) test_result2 = PropTest(control_sample, treatment_sample, 0.0, 1);

                // Stagewise test statistics
                stage1_test = test_result1.test_stat;
                stage2_test = test_result2.test_stat;

                pvalue2 = CombFunctionTestStat(stage1_test, stage2_test, weight, 1.0 - weight);

                outcome[3] = (pvalue2 <= alpha);        

            } 


        }
        // Normal and binary endpoints

        // Time-to-event endpoints
        if (endpoint_index == 3) {

            survival_data.os = fillvec(ntotal, 0.0);
            survival_data.os_local = fillvec(ntotal, 0.0);
            survival_data.os_local_censor = fillvec(ntotal, 0.0);
            survival_data.os_local_start = fillvec(ntotal, 0.0);
            survival_data.stratum = FillTreatmentIndicators(sample_size);

            // Patient enrollment time
            survival_data.start = Enrollment(ntotal, enrollment_period, enrollment_distribution, enrollment_parameter);    

            // Patient dropout time
            survival_data.dropout = Dropout(ntotal, 2, dropout_parameter); 

            // Populate the survival data frame and apply local transformations
            for (i = 0; i < ntotal; i++) {

                current_stratum = survival_data.stratum[i];

                survival_data.os[i] = Exponential(1, hazard_rates[current_stratum])[0];    

                survival_data.os_local_censor[i] = 0.0;

                // Local transformation for OS
                survival_data.os_local[i] = survival_data.os[i];
                if (survival_data.os[i] > survival_data.dropout[i]) {
                    survival_data.os_local[i] = survival_data.dropout[i];
                    survival_data.os_local_censor[i] = 1.0;
                }
                survival_data.os_local_start[i] = survival_data.os_local[i] + survival_data.start[i];
                if (survival_data.os_local_censor[i] == 1) survival_data.os_local_start[i] = -1.0;

            }

            // Outcomes
            // 0: Futility at IA1
            // 1: Significance at IA2
            // 2: Increase the number of events
            // 3: Significance at FA

            // Traditional analysis

            stratum_list.clear();
            stratum_list.push_back(0);        
            stratum_list.push_back(1);    
            look_time[0] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_fa1);

            // Control
            stratum_list.clear();
            stratum_list.push_back(0);
            outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[0]);

            // Treatment
            stratum_list.clear();
            stratum_list.push_back(1);
            outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[0]);

            test_result1 = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

            pvalue1 = test_result1.pvalue;

            trad_outcome = (pvalue1 <= alpha);        

            // Interim analysis 1

            stratum_list.clear();
            stratum_list.push_back(0);        
            stratum_list.push_back(1);    
            look_time[0] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_ia1);

            // Control
            stratum_list.clear();
            stratum_list.push_back(0);
            outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[0]);

            // Treatment
            stratum_list.clear();
            stratum_list.push_back(1);
            outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[0]);

            test_result1 = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);
            event_count[0] = EventCount(outcome_censor_control, outcome_censor_treatment);
            hr[0] = HazardRatio(outcome_censor_control, outcome_censor_treatment); 

            // Conditional power
            cp[0] = CondPower(test_result1.test_stat, event_count_ia1, event_count_fa1, 1, 1, 0.0, alpha);

            outcome[0] = (cp[0] <= futility_threshold);

            // Interim analysis 2

            if (outcome[0] < 1.0) {

                stratum_list.clear();
                stratum_list.push_back(0);        
                stratum_list.push_back(1);    
                look_time[1] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_ia2);

                // Control
                stratum_list.clear();
                stratum_list.push_back(0);
                outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[1]);

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(1);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[1]);

                test_result1 = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);
                event_count[1] = EventCount(outcome_censor_control, outcome_censor_treatment);
                hr[1] = HazardRatio(outcome_censor_control, outcome_censor_treatment); 

                // Conditional power
                cp[1] = CondPower(test_result1.test_stat, event_count_ia2, event_count_fa1, 1, 1, 0.0, alpha);

                pvalue1 = test_result1.pvalue;

            }

            if (outcome[0] < 1.0 && outcome[1] < 1.0) outcome[2] = (cp[1] > underpowered_zone[0] && cp[1] <= underpowered_zone[1]); 

            // Final analysis (before event count adjustment)
            if (outcome[0] < 1.0 && outcome[1] < 1.0 && outcome[2] < 1.0) {

                stratum_list.clear();
                stratum_list.push_back(0);        
                stratum_list.push_back(1);    
                look_time[2] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_fa1);

                // Control
                stratum_list.clear();
                stratum_list.push_back(0);
                outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(1);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

                test_result2 = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);
                event_count[2] = EventCount(outcome_censor_control, outcome_censor_treatment);
                hr[2] = HazardRatio(outcome_censor_control, outcome_censor_treatment); 

                pvalue2 = test_result2.pvalue;

                outcome[3] = (pvalue2 <= alpha);        

            } 

            // Final analysis (after event count adjustment)
            if (outcome[0] < 1.0 && outcome[1] < 1.0 && outcome[2] > 0.0) {

                // Updated event count
                sample_size_increase = UpdatedEventCount(test_result1.test_stat, event_count_ia2, event_count_fa1, event_count_fa2, 1, 1, 0.0, target_power, alpha); 

                stratum_list.clear();
                stratum_list.push_back(0);        
                stratum_list.push_back(1);    
                look_time[2] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_ia2 + sample_size_increase);

                // Control
                stratum_list.clear();
                stratum_list.push_back(0);
                outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(1);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

                test_result2 = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);
                event_count[2] = EventCount(outcome_censor_control, outcome_censor_treatment);
                hr[2] = HazardRatio(outcome_censor_control, outcome_censor_treatment); 

                // Stagewise test statistics
                stage1_test = test_result1.test_stat;
                stage2_test = (test_result2.test_stat * sqrt(event_count[2]) - test_result1.test_stat * sqrt(event_count[1])) / sqrt(event_count[2] - event_count[1]);

                pvalue2 = CombFunctionTestStat(stage1_test, stage2_test, weight, 1.0 - weight);

                outcome[3] = (pvalue2 <= alpha);        

            } 

        }
        // Time-to-event endpoints

        // Summary statistics for each simulation run

        sim_results(sim, 0) = outcome[0];
        sim_results(sim, 1) = outcome[1];
        sim_results(sim, 2) = outcome[2];
        sim_results(sim, 3) = outcome[3];

        sim_results(sim, 4) = look_time[0];
        sim_results(sim, 5) = look_time[1];
        sim_results(sim, 6) = look_time[2];

        sim_results(sim, 7) = cp[0];
        sim_results(sim, 8) = cp[1];
        sim_results(sim, 9) = event_count[0];
        sim_results(sim, 10) = event_count[1];
        sim_results(sim, 11) = event_count[2];
        sim_results(sim, 12) = sample_size_increase;
        sim_results(sim, 13) = hr[0];
        sim_results(sim, 14) = hr[1];
        sim_results(sim, 15) = hr[2];
        sim_results(sim, 16) = pvalue1;
        sim_results(sim, 17) = pvalue2;
        sim_results(sim, 18) = trad_outcome;


    }
    // End of simulations

    /*******************************************************************/

    // Return simulation results

    return List::create(Named("sim_results") = sim_results);

}
// End of ADSSModC


// [[Rcpp::export]]
List ADTreatSelC(const List &parameters_arg) {

    List parameters(parameters_arg);

    vector<int> sample_size = as<vector<int>>(parameters["sample_size"]);
    vector<int> sample_size_ia1 = as<vector<int>>(parameters["sample_size_ia1"]);
    vector<int> sample_size_ia2 = as<vector<int>>(parameters["sample_size_ia2"]);
    vector<int> sample_size_fa = as<vector<int>>(parameters["sample_size_fa"]);
    int max_sample_size = as<int>(parameters["max_sample_size"]);

    int event_count_ia1 = as<int>(parameters["event_count_ia1"]);
    int event_count_ia2 = as<int>(parameters["event_count_ia2"]);
    int event_count_fa = as<int>(parameters["event_count_fa"]);

    vector<double> dropout_parameter = as<vector<double>>(parameters["dropout_parameter"]);

    int endpoint_index = as<int>(parameters["endpoint_index"]);

    vector<double> means = as<vector<double>>(parameters["means"]);
    vector<double> sds = as<vector<double>>(parameters["sds"]);
    vector<double> rates = as<vector<double>>(parameters["rates"]);
    vector<double> hazard_rates = as<vector<double>>(parameters["hazard_rates"]);

    double futility_threshold = as<double>(parameters["futility_threshold"]);

    double enrollment_period = as<double>(parameters["enrollment_period"]);
    int enrollment_distribution = as<int>(parameters["enrollment_distribution"]);
    double enrollment_parameter = as<double>(parameters["enrollment_parameter"]);

    int nsims = as<int>(parameters["nsims"]);
    double alpha = as<double>(parameters["alpha"]);

    /*******************************************************************/

    double pvalue, ad_outcome, max_cp;

    int i, j, sim, current_stratum, ntotal, narms = sample_size.size(), select_flag;

    ntotal = SumVecInt(sample_size);

    AllSurvivalData survival_data;

    OutcomeCensor outcome_censor, outcome_censor_control, outcome_censor_treatment;

    TestResult test_result;

    vector<double> look_time(3), event_count(3), dropout, hr(narms), cp1(narms - 1), cp2(narms - 1), overall_control_sample, control_sample, treatment_sample, temp_vec, pvalue_trad(narms - 1), trad_outcome(narms - 1);

    vector<int> stratum_list, sample_list, futility_flag(narms - 1);

    NumericMatrix sim_results(nsims, 10 + narms), overall_treatment_sample(narms - 1, max_sample_size); 

    /*******************************************************************/

    for (sim = 0; sim < nsims; sim++) {

        hr = fillvec(narms, 0.0);
        cp1 = fillvec(narms - 1, 0.0);
        cp2 = fillvec(narms - 1, 0.0);
        trad_outcome = fillvec(narms - 1, 0.0);

        look_time = fillvec(3, 0.0);
        event_count = fillvec(3, 0.0);
        pvalue_trad = fillvec(narms - 1, 0.0);

        pvalue = 0.0;

        select_flag = -1;
        max_cp = 0.0;

        ad_outcome = 0.0;

        // Normal and binary endpoints
        if (endpoint_index == 1 || endpoint_index == 2) {

            // Control samples
            if (endpoint_index == 1) overall_control_sample = Normal(sample_size_fa[0], means[0], sds[0]);
            if (endpoint_index == 2) overall_control_sample = Binary(sample_size_fa[0], rates[0]);

            // Treatment samples
            for (i = 0; i < narms - 1; i++) {

                if (endpoint_index == 1) temp_vec = Normal(sample_size_fa[i + 1], means[i + 1], sds[i + 1]);
                if (endpoint_index == 2) temp_vec = Binary(sample_size_fa[i + 1], rates[i + 1]);

                for (j = 0; j < sample_size_fa[i + 1]; j++) overall_treatment_sample(i, j) = temp_vec[j];

            }

            // Interim analysis 1      

            control_sample = ExtractSamples(overall_control_sample, 0, sample_size_ia1[0]);

            for (i = 0; i < narms - 1; i++) {

                treatment_sample.clear();
                for (j = 0; j < sample_size_ia1[i + 1]; j++) treatment_sample.push_back(overall_treatment_sample(i, j));

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, 1);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, 1);

                // Conditional power
                cp1[i] = CondPower(test_result.test_stat, sample_size_ia1[0] + sample_size_ia1[i + 1], sample_size_fa[0] + sample_size_fa[i + 1], 1, 1, 0.0, alpha);

                // Futility stopping rule
                futility_flag[i] = (cp1[i] <= futility_threshold);

            }

// printVectorDouble("cp1", cp1);
// printVectorInt("futility_flag", futility_flag);

            // Interim analysis 2

            control_sample = ExtractSamples(overall_control_sample, 0, sample_size_ia2[0]);

            for (i = 0; i < narms - 1; i++) {

                treatment_sample.clear();
                for (j = 0; j < sample_size_ia2[i + 1]; j++) treatment_sample.push_back(overall_treatment_sample(i, j));

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, 1);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, 1);

                // Conditional power
                cp2[i] = CondPower(test_result.test_stat, sample_size_ia2[0] + sample_size_ia2[i + 1], sample_size_fa[0] + sample_size_fa[i + 1], 1, 1, 0.0, alpha);

                // Treatment selection rule
                if (futility_flag[i] < 1.0 && cp2[i] > max_cp) {
                    max_cp = cp2[i];
                    select_flag = i;    
                }

            }

// printVectorDouble("cp2", cp2);
// cerr<<"select_flag: "<<select_flag<<endl;

            // Final analysis

            if (select_flag > -1) {

                control_sample = ExtractSamples(overall_control_sample, 0, sample_size_fa[0]);

                treatment_sample.clear();
                for (j = 0; j < sample_size_fa[select_flag + 1]; j++) treatment_sample.push_back(overall_treatment_sample(select_flag, j));

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, 1);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, 1);

                ad_outcome = (test_result.pvalue <= alpha / (narms - 1.0));        

            }

            // Traditional analysis

            control_sample = ExtractSamples(overall_control_sample, 0, sample_size_fa[0]);

            for (i = 0; i < narms - 1; i++) {

                treatment_sample.clear();
                for (j = 0; j < sample_size_fa[i + 1]; j++) treatment_sample.push_back(overall_treatment_sample(i, j));

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, 1);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, 1);

                pvalue_trad[i] = test_result.pvalue;

                trad_outcome[i] = (pvalue_trad[i] <= alpha && futility_flag[i] < 1.0);        

            }



        }
        // Normal and binary endpoints

        // Time-to-event endpoints
        if (endpoint_index == 3) {

            survival_data.os = fillvec(ntotal, 0.0);
            survival_data.os_local = fillvec(ntotal, 0.0);
            survival_data.os_local_censor = fillvec(ntotal, 0.0);
            survival_data.os_local_start = fillvec(ntotal, 0.0);
            survival_data.stratum = FillTreatmentIndicators(sample_size);

            // Patient enrollment time
            survival_data.start = Enrollment(ntotal, enrollment_period, enrollment_distribution, enrollment_parameter);    

            // Patient dropout time
            survival_data.dropout = Dropout(ntotal, 2, dropout_parameter); 

            // Populate the survival data frame and apply local transformations
            for (i = 0; i < ntotal; i++) {

                current_stratum = survival_data.stratum[i];

                survival_data.os[i] = Exponential(1, hazard_rates[current_stratum])[0];    

                survival_data.os_local_censor[i] = 0.0;

                // Local transformation for OS
                survival_data.os_local[i] = survival_data.os[i];
                if (survival_data.os[i] > survival_data.dropout[i]) {
                    survival_data.os_local[i] = survival_data.dropout[i];
                    survival_data.os_local_censor[i] = 1.0;
                }
                survival_data.os_local_start[i] = survival_data.os_local[i] + survival_data.start[i];
                if (survival_data.os_local_censor[i] == 1) survival_data.os_local_start[i] = -1.0;

            }

            // Interim analysis 1

            stratum_list.clear();
            for (i = 0; i < narms; i++) stratum_list.push_back(i);        
            look_time[0] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_ia1);

            // Control
            stratum_list.clear();
            stratum_list.push_back(0);
            outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[0]);

            for (i = 0; i < narms - 1; i++) {

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(i + 1);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[0]);
                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

                // Conditional power
                cp1[i] = CondPower(test_result.test_stat, event_count_ia1, event_count_fa, 1, 1, 0.0, alpha);

                // Futility stopping rule
                futility_flag[i] = (cp1[i] <= futility_threshold);

            }

            // Interim analysis 1

            stratum_list.clear();
            for (i = 0; i < narms; i++) stratum_list.push_back(i);        
            look_time[1] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_ia2);

            // Control
            stratum_list.clear();
            stratum_list.push_back(0);
            outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[1]);

            for (i = 0; i < narms - 1; i++) {

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(i + 1);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[1]);
                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

                // Conditional power
                cp2[i] = CondPower(test_result.test_stat, event_count_ia2, event_count_fa, 1, 1, 0.0, alpha);

                // Treatment selection rule
                if (futility_flag[i] < 1.0 && cp2[i] > max_cp) {
                    max_cp = cp2[i];
                    select_flag = i;    
                }

            }

            // Final analysis

            stratum_list.clear();
            for (i = 0; i < narms; i++) stratum_list.push_back(i);        
            look_time[2] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_fa);

            if (select_flag > -1) {

                // Control
                stratum_list.clear();
                stratum_list.push_back(0);
                outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(select_flag + 1);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);
    
                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

                ad_outcome = (test_result.pvalue <= alpha / (narms - 1.0));        

            }

            // Traditional analysis

            // Control
            stratum_list.clear();
            stratum_list.push_back(0);
            outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

            for (i = 0; i < narms - 1; i++) {

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(i + 1);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);
    
                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

                pvalue_trad[i] = test_result.pvalue;

                trad_outcome[i] = (pvalue_trad[i] <= alpha && futility_flag[i] < 1.0);        

            }

        }
        // Time-to-event endpoints

        // Summary statistics for each simulation run

        sim_results(sim, 0) = ad_outcome;
        sim_results(sim, 1) = select_flag;
        for (i = 0; i < narms - 1; i++) sim_results(sim, 2 + i) = trad_outcome[i];
        for (i = 0; i < narms - 1; i++) sim_results(sim, 2 + narms - 1 + i) = futility_flag[i];

    }
    // End of simulations

    /*******************************************************************/

    // Return simulation results

    return List::create(Named("sim_results") = sim_results);

}
// End of ADTreatSelC

// [[Rcpp::export]]
List ADPopSelC(const List &parameters_arg) {

    List parameters(parameters_arg);

    vector<int> sample_size_pop = as<vector<int>>(parameters["sample_size_pop"]);
    vector<int> sample_size = as<vector<int>>(parameters["sample_size"]);
    vector<int> sample_size_ia1 = as<vector<int>>(parameters["sample_size_ia1"]);
    vector<int> sample_size_ia2 = as<vector<int>>(parameters["sample_size_ia2"]);
    vector<int> sample_size_fa = as<vector<int>>(parameters["sample_size_fa"]);

    int event_count_ia1 = as<int>(parameters["event_count_ia1"]);
    int event_count_ia2 = as<int>(parameters["event_count_ia2"]);
    vector<int> event_count_fa = as<vector<int>>(parameters["event_count_fa"]);

    vector<double> dropout_parameter = as<vector<double>>(parameters["dropout_parameter"]);

    int endpoint_index = as<int>(parameters["endpoint_index"]);
    int max_sample_size = as<int>(parameters["max_sample_size"]);

    vector<double> means = as<vector<double>>(parameters["means"]);
    vector<double> sds = as<vector<double>>(parameters["sds"]);
    vector<double> rates = as<vector<double>>(parameters["rates"]);
    vector<double> hazard_rates = as<vector<double>>(parameters["hazard_rates"]);

    double futility_threshold = as<double>(parameters["futility_threshold"]);
    double influence = as<double>(parameters["influence"]);
    double interaction = as<double>(parameters["interaction"]);

    double enrollment_period = as<double>(parameters["enrollment_period"]);
    int enrollment_distribution = as<int>(parameters["enrollment_distribution"]);
    double enrollment_parameter = as<double>(parameters["enrollment_parameter"]);

    int nsims = as<int>(parameters["nsims"]);
    double alpha = as<double>(parameters["alpha"]);

    /*******************************************************************/

    double cp, futility_flag, effect_size_minus, effect_size_plus, hr_minus, hr_plus, temp;

    int i, j, sim, current_stratum, ntotal = sample_size[0] + sample_size[1], narms = 4;

    AllSurvivalData survival_data;

    OutcomeCensor outcome_censor, outcome_censor_control, outcome_censor_treatment;

    TestResult test_result;

    vector<double> look_time(4), event_count(4), dropout, overall_control_sample, overall_treatment_sample, control_sample, treatment_sample, temp_vec, trad_outcome(2), hypothesis_selection_outcome(3), ad_outcome(2), sample0(max_sample_size), sample1(max_sample_size), sample2(max_sample_size), sample3(max_sample_size), pvalue(2), outcome(2);

    vector<int> stratum_list, sample_list;

    NumericMatrix sim_results(nsims, 19); 

    /*******************************************************************/

    for (sim = 0; sim < nsims; sim++) {

        ad_outcome = fillvec(2, 0.0);
        trad_outcome = fillvec(2, 0.0);
        look_time = fillvec(4, 0.0);
        event_count = fillvec(4, 0.0);
        hypothesis_selection_outcome = fillvec(3, 0.0);

        cp = 0.0;

        // Normal and binary endpoints
        if (endpoint_index == 1 || endpoint_index == 2) {

            // Overall samples
            for (j = 0; j < 4; j++) {

                if (endpoint_index == 1) temp_vec = Normal(sample_size_fa[j], means[j], sds[j]);
                if (endpoint_index == 2) temp_vec = Binary(sample_size_fa[j], rates[j]);
                
                for (i = 0; i < sample_size_fa[j]; i++) {
                    if (j == 0) sample0[i] = temp_vec[i];
                    if (j == 1) sample1[i] = temp_vec[i];
                    if (j == 2) sample2[i] = temp_vec[i];
                    if (j == 3) sample3[i] = temp_vec[i];
                }

            }

            // Interim analysis 1      

            control_sample = CombineVec(ExtractSamples(sample0, 0, sample_size_ia1[0]), ExtractSamples(sample1, 0, sample_size_ia1[1]));
            treatment_sample = CombineVec(ExtractSamples(sample2, 0, sample_size_ia1[2]), ExtractSamples(sample3, 0, sample_size_ia1[3]));

            if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, 1);
            if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, 1);

            // Conditional power
            cp = CondPower(test_result.test_stat, SumVecInt(sample_size_ia1), SumVecInt(sample_size_fa), 1, 1, 0.0, alpha);

            // Futility stopping rule
            futility_flag = (cp <= futility_threshold);

            // Interim analysis 2

            // Biomarker-negative subpopulation

            control_sample = ExtractSamples(sample0, 0, sample_size_ia2[0]);
            treatment_sample = ExtractSamples(sample2, 0, sample_size_ia2[2]);

            effect_size_minus = ComputeEffectSize(control_sample, treatment_sample, endpoint_index, 1);

            // Biomarker-positive subpopulation

            control_sample = ExtractSamples(sample1, 0, sample_size_ia2[1]);
            treatment_sample = ExtractSamples(sample3, 0, sample_size_ia2[3]);

            effect_size_plus = ComputeEffectSize(control_sample, treatment_sample, endpoint_index, 1);

            hypothesis_selection_outcome = HypothesisSelection(effect_size_minus, effect_size_plus, influence, interaction); 

            // Final analysis

            // All patients are included in the final analysis
            if (hypothesis_selection_outcome[0] > 0.0 || hypothesis_selection_outcome[2] > 0.0) {

                // Overall population

                control_sample = CombineVec(ExtractSamples(sample0, 0, sample_size_fa[0]), ExtractSamples(sample1, 0, sample_size_fa[1]));
                treatment_sample = CombineVec(ExtractSamples(sample2, 0, sample_size_fa[2]), ExtractSamples(sample3, 0, sample_size_fa[3]));

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, 1);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, 1);

                ad_outcome[0] = (test_result.pvalue <= alpha && futility_flag < 1.0);

                // Biomarker-positive subpopulation

                control_sample = ExtractSamples(sample1, 0, sample_size_fa[1]);
                treatment_sample = ExtractSamples(sample3, 0, sample_size_fa[3]);

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, 1);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, 1);

                ad_outcome[1] = (test_result.pvalue <= alpha && futility_flag < 1.0);

            } 

            // Only biomarker-positive patients are included in the final analysis
            if (hypothesis_selection_outcome[1] > 0.0) {

                control_sample = ExtractSamples(sample1, 0, sample_size_fa[1]);
                treatment_sample = ExtractSamples(sample3, 0, sample_size_fa[3]);

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, 1);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, 1);

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

                ad_outcome[1] = (test_result.pvalue <= alpha && futility_flag < 1.0);

            } 

            // Traditional analysis

            // Overall population

            control_sample = CombineVec(ExtractSamples(sample0, 0, sample_size_fa[0]), ExtractSamples(sample1, 0, sample_size_fa[1]));
            treatment_sample = CombineVec(ExtractSamples(sample2, 0, sample_size_fa[2]), ExtractSamples(sample3, 0, sample_size_fa[3]));

            if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, 1);
            if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, 1);

            trad_outcome[0] = (test_result.pvalue <= alpha && futility_flag < 1.0);        

            // Biomarker-positive subpopulation

            control_sample = ExtractSamples(sample1, 0, sample_size_fa[1]);
            treatment_sample = ExtractSamples(sample3, 0, sample_size_fa[3]);

            if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, 1);
            if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, 1);

            trad_outcome[1] = (test_result.pvalue <= alpha && futility_flag < 1.0);        


        }
        // Normal and binary endpoints

        // Time-to-event endpoints
        if (endpoint_index == 3) {

            survival_data.os = fillvec(ntotal, 0.0);
            survival_data.os_local = fillvec(ntotal, 0.0);
            survival_data.os_local_censor = fillvec(ntotal, 0.0);
            survival_data.os_local_start = fillvec(ntotal, 0.0);
            survival_data.stratum = FillTreatmentIndicators(sample_size_pop);

            // Patient enrollment time
            survival_data.start = Enrollment(ntotal, enrollment_period, enrollment_distribution, enrollment_parameter);    

            // Patient dropout time
            survival_data.dropout = Dropout(ntotal, 2, dropout_parameter); 

            // Populate the survival data frame and apply local transformations
            for (i = 0; i < ntotal; i++) {

                current_stratum = survival_data.stratum[i];

                survival_data.os[i] = Exponential(1, hazard_rates[current_stratum])[0];    

                survival_data.os_local_censor[i] = 0.0;

                // Local transformation for OS
                survival_data.os_local[i] = survival_data.os[i];
                if (survival_data.os[i] > survival_data.dropout[i]) {
                    survival_data.os_local[i] = survival_data.dropout[i];
                    survival_data.os_local_censor[i] = 1.0;
                }
                survival_data.os_local_start[i] = survival_data.os_local[i] + survival_data.start[i];
                if (survival_data.os_local_censor[i] == 1) survival_data.os_local_start[i] = -1.0;

            }

            // Interim analysis 1

            stratum_list.clear();
            for (i = 0; i < narms; i++) stratum_list.push_back(i);        
            look_time[0] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_ia1);

            // Control
            stratum_list.clear();
            stratum_list.push_back(0);
            stratum_list.push_back(1);
            outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[0]);

            // Treatment
            stratum_list.clear();
            stratum_list.push_back(2);
            stratum_list.push_back(3);
            outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[0]);
            test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);
            event_count[0] = EventCount(outcome_censor_control, outcome_censor_treatment);

            // Conditional power
            cp = CondPower(test_result.test_stat, event_count_ia1, event_count_fa[0], 1, 1, 0.0, alpha);

            // Futility stopping rule
            futility_flag = (cp <= futility_threshold);

            // Interim analysis 2

            stratum_list.clear();
            for (i = 0; i < narms; i++) stratum_list.push_back(i);        
            look_time[1] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_ia2);

            // Biomarker-negative subpopulation

            // Control
            stratum_list.clear();
            stratum_list.push_back(0);
            outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[1]);

            // Treatment
            stratum_list.clear();
            stratum_list.push_back(2);
            outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[1]);

            hr_minus = HazardRatio(outcome_censor_control, outcome_censor_treatment);

            // Biomarker-positive subpopulation

            // Control
            stratum_list.clear();
            stratum_list.push_back(1);
            outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[1]);

            // Treatment
            stratum_list.clear();
            stratum_list.push_back(3);
            outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[1]);
            event_count[1] = EventCount(outcome_censor_control, outcome_censor_treatment);

            hr_plus = HazardRatio(outcome_censor_control, outcome_censor_treatment);

            effect_size_minus = -log(hr_minus);
            effect_size_plus = -log(hr_plus);

            hypothesis_selection_outcome = HypothesisSelection(effect_size_minus, effect_size_plus, influence, interaction); 

            // Final analysis

            // Only in the overall population

            if (hypothesis_selection_outcome[0] > 0.0) {

                stratum_list.clear();
                for (i = 0; i < narms; i++) stratum_list.push_back(i);        
                look_time[2] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_fa[0]);

                // Overall population

                // Control
                stratum_list.clear();
                stratum_list.push_back(0);
                stratum_list.push_back(1);
                outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(2);
                stratum_list.push_back(3);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

                ad_outcome[0] = (test_result.pvalue <= alpha / 2.0 && futility_flag < 1.0);

            } 

            // Only in the biomarker-positive population
            if (hypothesis_selection_outcome[1] > 0.0) {

                stratum_list.clear();
                stratum_list.push_back(1);
                stratum_list.push_back(3);
                look_time[3] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_fa[1]);

                // Control
                stratum_list.clear();
                stratum_list.push_back(1);
                outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[3]);

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(3);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[3]);
                event_count[3] = EventCount(outcome_censor_control, outcome_censor_treatment);

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

                ad_outcome[1] = (test_result.pvalue <= alpha / 2.0 && futility_flag < 1.0);

            } 

            // In both populations

            if (hypothesis_selection_outcome[2] > 0.0) {

                stratum_list.clear();
                for (i = 0; i < narms; i++) stratum_list.push_back(i);        
                look_time[2] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_fa[0]);

                // Overall population

                // Control
                stratum_list.clear();
                stratum_list.push_back(0);
                stratum_list.push_back(1);
                outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(2);
                stratum_list.push_back(3);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

                pvalue[0] = test_result.pvalue;

                // Biomarker-positive subpopulation

                // Control
                stratum_list.clear();
                stratum_list.push_back(1);
                outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(3);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);
                event_count[2] = EventCount(outcome_censor_control, outcome_censor_treatment);

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

                pvalue[1] = test_result.pvalue;

                outcome = HochbergOutcome(pvalue, alpha);

                ad_outcome[0] = (outcome[0] > 0.0 && futility_flag < 1.0);
                ad_outcome[1] = (outcome[1] > 0.0 && futility_flag < 1.0);

            } 


            // Traditional analysis

            // Overall population

            stratum_list.clear();
            for (i = 0; i < narms; i++) stratum_list.push_back(i);        
            temp = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_fa[0]);

            // Control
            stratum_list.clear();
            stratum_list.push_back(0);
            stratum_list.push_back(1);
            outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, temp);

            // Treatment
            stratum_list.clear();
            stratum_list.push_back(2);
            stratum_list.push_back(3);
            outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, temp);

            test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

            trad_outcome[0] = (test_result.pvalue <= alpha && futility_flag < 1.0);        

            // Biomarker-positive subpopulation

            stratum_list.clear();
            stratum_list.push_back(1);
            stratum_list.push_back(3);
            temp = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_fa[1]);

            // Control
            stratum_list.clear();
            stratum_list.push_back(1);
            outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, temp);

            // Treatment
            stratum_list.clear();
            stratum_list.push_back(3);
            outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, temp);

            test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

            trad_outcome[1] = (test_result.pvalue <= alpha && futility_flag < 1.0); 

        }
        // Time-to-event endpoints

        // Summary statistics for each simulation run

        sim_results(sim, 0) = trad_outcome[0];
        sim_results(sim, 1) = trad_outcome[1];
        sim_results(sim, 2) = ad_outcome[0];
        sim_results(sim, 3) = ad_outcome[1];

        sim_results(sim, 4) = cp;
        sim_results(sim, 5) = futility_flag;
        sim_results(sim, 6) = effect_size_minus;
        sim_results(sim, 7) = effect_size_plus;
        sim_results(sim, 8) = hypothesis_selection_outcome[0];
        sim_results(sim, 9) = hypothesis_selection_outcome[1];
        sim_results(sim, 10) = hypothesis_selection_outcome[2];

        sim_results(sim, 11) = look_time[0];
        sim_results(sim, 12) = look_time[1];
        sim_results(sim, 13) = look_time[2];
        sim_results(sim, 14) = look_time[3];

        sim_results(sim, 15) = event_count[1];
        sim_results(sim, 16) = event_count[2];

    }
    // End of simulations

    /*******************************************************************/

    // Return simulation results

    return List::create(Named("sim_results") = sim_results);

}
// End of ADPopSelC

// [[Rcpp::export]]
List FutRuleC(const List &parameters_arg) {

    List parameters(parameters_arg);

    vector<int> sample_size = as<vector<int>>(parameters["sample_size"]);
    vector<int> interim = as<vector<int>>(parameters["interim"]);
    vector<int> final = as<vector<int>>(parameters["final"]);
    vector<double> dropout_parameter = as<vector<double>>(parameters["dropout_parameter"]);

    vector<double> means = as<vector<double>>(parameters["means"]);
    vector<double> sds = as<vector<double>>(parameters["sds"]);
    vector<double> rates = as<vector<double>>(parameters["rates"]);
    vector<double> hazard_rates = as<vector<double>>(parameters["hazard_rates"]);

    int endpoint_index = as<int>(parameters["endpoint_index"]);

    int nsims = as<int>(parameters["nsims"]);
    double alpha = as<double>(parameters["alpha"]);

    double enrollment_period = as<double>(parameters["enrollment_period"]);
    int enrollment_distribution = as<int>(parameters["enrollment_distribution"]);
    double enrollment_parameter = as<double>(parameters["enrollment_parameter"]);

    /*******************************************************************/

    int i, sim, narms = sample_size.size(), n, current_stratum;

    double look_time;

    n = SumVecInt(sample_size);

    AllSurvivalData survival_data;

    OutcomeCensor outcome_censor, outcome_censor_control, outcome_censor_treatment;

    TestResult test_result;

    vector<double> control_sample, treatment_sample, cp(narms - 1);
    vector<int> stratum_list;

    NumericMatrix sim_results(nsims, narms - 1);

    /*******************************************************************/

    for (sim = 0; sim < nsims; sim++) {

        // Normal and binary endpoints
        if (endpoint_index == 1 || endpoint_index == 2) {

            // Control sample
            if (endpoint_index == 1) control_sample = Normal(interim[0], means[0], sds[0]);
            if (endpoint_index == 2) control_sample = Binary(interim[0], rates[0]);

            // Compute conditional power

            for (i = 0; i < narms - 1; i++) {

                // Normal endpoint
                if (endpoint_index == 1) {
                    treatment_sample = Normal(interim[i + 1], means[i + 1], sds[i + 1]);
                    test_result = TTest(control_sample, treatment_sample, 0.0, 1);
                }

                // Binary endpoint
                if (endpoint_index == 2) {
                    treatment_sample = Binary(interim[i + 1], rates[i + 1]);
                    test_result = PropTest(control_sample, treatment_sample, 0.0, 1);
                }

                // Conditional power
                cp[i] = CondPower(test_result.test_stat, interim[0], final[0], 1, 1, 0.0, alpha);

                sim_results(sim, i) = cp[i];

            }

        }

        // Time-to-event endpoints
        if (endpoint_index == 3) {

            survival_data.os = fillvec(n, 0.0);
            survival_data.os_local = fillvec(n, 0.0);
            survival_data.os_local_censor = fillvec(n, 0.0);
            survival_data.os_local_start = fillvec(n, 0.0);
            survival_data.stratum = FillTreatmentIndicators(sample_size);

            // Patient enrollment time
            survival_data.start = Enrollment(n, enrollment_period, enrollment_distribution, enrollment_parameter);    

            // Patient dropout time
            survival_data.dropout = Dropout(n, 2, dropout_parameter); 

            // Populate the survival data frame and apply local transformations
            for (i = 0; i < n; i++) {

                current_stratum = survival_data.stratum[i];

                survival_data.os[i] = Exponential(1, hazard_rates[current_stratum])[0];    

                survival_data.os_local_censor[i] = 0.0;

                // Local transformation for OS
                survival_data.os_local[i] = survival_data.os[i];
                if (survival_data.os[i] > survival_data.dropout[i]) {
                    survival_data.os_local[i] = survival_data.dropout[i];
                    survival_data.os_local_censor[i] = 1.0;
                }
                survival_data.os_local_start[i] = survival_data.os_local[i] + survival_data.start[i];
                if (survival_data.os_local_censor[i] == 1) survival_data.os_local_start[i] = -1.0;

            }

            // Interim analysis

            stratum_list.clear();
            for (i = 0; i < narms; i++) stratum_list.push_back(i);        
            look_time = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, interim[0]);

            // Compute conditional power

            for (i = 0; i < narms - 1; i++) {

                // Control
                stratum_list.clear();
                stratum_list.push_back(0);
                outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time);

                // Treatment
                stratum_list.clear();
                stratum_list.push_back(i + 1);
                outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time);

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, 1);

                // Conditional power
                cp[i] = CondPower(test_result.test_stat, interim[0], final[0], 1, 1, 0.0, alpha);

                sim_results(sim, i) = cp[i];

            }


        }
        // Time-to-event endpoints



    }
    // End of simulations

    /*******************************************************************/

    // Return simulation results

    return List::create(Named("sim_results") = sim_results);

}
// End of FutRuleC

// [[Rcpp::export]]
List EventPredEventCount(const List &parameters_arg) {

    List parameters(parameters_arg);

    vector<double> event = as<vector<double>>(parameters["event"]);
    vector<double> dropout = as<vector<double>>(parameters["dropout"]);
    vector<double> enrollment = as<vector<double>>(parameters["enrollment"]);
    vector<double> time = as<vector<double>>(parameters["time"]);
    vector<double> event_sample = as<vector<double>>(parameters["event_sample"]);
    vector<double> dropout_sample = as<vector<double>>(parameters["dropout_sample"]);
    vector<double> enrollment_sample = as<vector<double>>(parameters["enrollment_sample"]);
    double time_point = as<double>(parameters["time_point"]);
    double interim_analysis = as<double>(parameters["interim_analysis"]);

    /*******************************************************************/

    int i, nsims = event_sample.size(), n = event.size(), sim;

    vector<double> sim_results(nsims), complete_event;

    // Length of maximum patient follow up
    double rel_event_time, rel_dropout_time, current_enrollment_time, event_indicator;

    /*******************************************************************/

    // Create a complete data set using posterior parameters and compute the number of events at the current time point

    for (sim = 0; sim < nsims; sim++) {      

       complete_event.clear(); 

       // Patients enrolled prior to the interim analysis  
       for (i = 0; i < n; i++) {

            if (event[i] > 0.0) complete_event.push_back(1.0);

            if (event[i] < 1.0 && dropout[i] < 1.0) {

              // Relative event time 
              rel_event_time = time[i] + Exponential(1, event_sample[sim])[0];

              // Relative dropout time 
              rel_dropout_time = time[i] + Exponential(1, dropout_sample[sim])[0];

              event_indicator = (rel_dropout_time > rel_event_time && rel_event_time < time_point - enrollment[i]);

              complete_event.push_back(event_indicator);


            } 

       }

       // Patients enrolled after the interim analysis  

       current_enrollment_time = interim_analysis;

       do {

            current_enrollment_time += Exponential(1, enrollment_sample[sim])[0];

            // Relative event time 
            rel_event_time = Exponential(1, event_sample[sim])[0];

            // Relative dropout time 
            rel_dropout_time = Exponential(1, dropout_sample[sim])[0];

            event_indicator = (rel_dropout_time > rel_event_time && rel_event_time < time_point - current_enrollment_time);

            complete_event.push_back(event_indicator);


       }
       while (current_enrollment_time <= time_point);

       sim_results[sim] = SumVec(complete_event);

   }

    /*******************************************************************/

    // Return simulation results

    return List::create(Named("sim_results") = sim_results);

}
// End of EventPredEventCount



