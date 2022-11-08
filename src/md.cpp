#include <Rcpp.h>
#include <RcppNumerical.h>

using namespace Numer;
using namespace std;
using namespace Rcpp;
using namespace Eigen;

#include "medstruct.h"
#include "medsupport.h"
#include "meddist.h"
#include "medstattest.h"
#include "medmult.h"
#include "varcuboid.h"

#include <exception>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <math.h> 
#include <iomanip>
#include <fstream>
#include <string>

//TODO: Move to calculation structures for ADRand
int n_models = 4, direction_index_local = 1;
//TODO: Move to calculation structures for ADRand
double final_gradient;



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

            if (pvalue[0] > alpha && pvalue[1] <= alpha / 2) {      // # nocov start
                outcome[0] = 0;
                outcome[1] = 1;
            }                                                       // # nocov end

        }

    }

    return outcome;

}

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

double ComputeEffectSize(const std::vector<double> &sample1, const std::vector<double> &sample2, const int &endpoint_distribution, const int &endpoint_test, const int &direction) {

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
        if (direction == 2) effect_size = -effect_size;

    }

    // Binary distribution
    if (endpoint_distribution == 2) {

        double estimatex = ComputeRate(sample1);
        double estimatey = ComputeRate(sample2);
        double ave = (estimatex+estimatey)/2.0;

        // Without pooled variance
        if (endpoint_test == 1) effect_size = (estimatey - estimatex) / sqrt(estimatex*(1.0-estimatex)+ estimatey*(1.0-estimatey)); 

        // With pooled variance
        if (endpoint_test == 2) effect_size = (estimatey - estimatex) / sqrt(ave*(1.0-ave)); 

        if (direction == 2) effect_size = -effect_size;

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

        n = 0;      // # nocov

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

        // # nocov start
        // Both populations are always selected if the interaction condition is not specified
        outcome3 = 1.0;
        // # nocov end

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

// REM, 43
double Logit(const double &x) {
    return log(x /(1.0 - x));
}

// REM, 47
// Compute the mean of a vector
double MeanVec(const vector<double> &vec)
{
   double mean_vec = 0.0;

   if (vec.size() > 0) mean_vec = accumulate(vec.begin(), vec.end(), 0.0) / vec.size();

   return(mean_vec);
}

//REM, 146
double DoseResponseFunction(const double &x, const int &model, const vector<double> &beta, const double &direction_index) {

    double y = 0.0, den;

    // Linear model
    if (model == 1) {
        y = beta[0] + beta[1] * x;
    }

    // Exponential model
    if (model == 2) {
        y = beta[0] + beta[1] * (exp(x / beta[2]) - 1.0);
    }

    // Emax model
    if (model == 3) {
        y = beta[0] + beta[1] * x / (beta[2] + x);
    }

    // Logistic model
    if (model == 4) {
        den = 1.0 + exp((beta[2] - x) / beta[3]);
        y = beta[0] + beta[1] / den;
    }

    return y;

}

//REM, 176
// Compute the model parameters to match the placebo and maximum effects
vector<double> StandardDRFunction(const int &model_index, const double &placebo_effect, const double &max_effect, const double &max_dose, const vector<double> &parameters) {

    int i;
    vector<double> coef;    

    // Linear model
    if (model_index == 1) {
        vector<double> local_coef(2);
        local_coef[0] = placebo_effect;
        local_coef[1] = max_effect / max_dose;
        for (i = 0; i < 2; i++) coef.push_back(local_coef[i]);
    }

    // Exponential model
    if (model_index == 2) {
        vector<double> local_coef(3);
        local_coef[0] = placebo_effect;
        local_coef[1] = max_effect / (exp(max_dose / parameters[0]) - 1.0);
        local_coef[2] = parameters[0];
        for (i = 0; i < 3; i++) coef.push_back(local_coef[i]);
    }

    // Emax model
    if (model_index == 3) {
        vector<double> local_coef(3);
        local_coef[0] = placebo_effect;
        local_coef[1] = max_effect * (parameters[0] + max_dose) / max_dose;
        local_coef[2] = parameters[0];
        for (i = 0; i < 3; i++) coef.push_back(local_coef[i]);
    }

    // Logistic model
    if (model_index == 4) {
        vector<double> local_coef(4);
        vector<double> temp_local_coef(4);
        temp_local_coef[0] = 0.0;
        temp_local_coef[1] = 1.0;
        temp_local_coef[2] = parameters[0];
        temp_local_coef[3] = parameters[1];

        double temp =  max_effect / (DoseResponseFunction(max_dose, 4, temp_local_coef, 1.0) - DoseResponseFunction(0.0, 4, temp_local_coef, 1.0));
        local_coef[0] = placebo_effect - temp * DoseResponseFunction(0.0, 4, temp_local_coef, 1.0);
        local_coef[1] = temp;
        local_coef[2] = parameters[0];
        local_coef[3] = parameters[1];
        for (i = 0; i < 4; i++) coef.push_back(local_coef[i]);
    }

    return coef; 

}

//REM, 229
// Compute the model parameters to match the placebo and maximum effects
vector<double> ComputeDoseResponseFunctionParameters(const int &model, const double &placebo_effect, const double &max_effect, const double &max_dose, const vector<double> &nonlinear_parameters) {

    vector<double> coef(5), temp_coef(5);
    double temp = 0.0;

    // Exponential model
    if (model == 2) {
        coef[0] = placebo_effect;
        coef[1] = max_effect / (exp(max_dose / nonlinear_parameters[0]) - 1.0);
        coef[2] = nonlinear_parameters[0];
    }

    // Emax model
    if (model == 3) {
        coef[0] = placebo_effect;
        coef[1] = max_effect * (nonlinear_parameters[0] + max_dose) / max_dose;
        coef[2] = nonlinear_parameters[0];
    }

    // Logistic model
    if (model == 4) {
        temp_coef[0] = 0.0;
        temp_coef[1] = 1.0;
        temp_coef[2] = nonlinear_parameters[0];
        temp_coef[3] = nonlinear_parameters[1];
        temp =  max_effect / (DoseResponseFunction(max_dose, 4, temp_coef, 1) - DoseResponseFunction(0.0, 4, temp_coef, 1));
        coef[0] = placebo_effect - temp * DoseResponseFunction(0.0, 4, temp_coef, 1);
        coef[1] = temp;
        coef[2] = nonlinear_parameters[0];
        coef[3] = nonlinear_parameters[1];
    }

    return coef;

}

//REM, 269
NumericMatrix OptContrast(const NumericMatrix &user_specified, const vector<int> &model_index, const double &direction_index, const vector<double> &dose_levels, const vector<double> &var) {

    int i, j, n_doses = dose_levels.size(), n_models = model_index.size();
    double max_dose = dose_levels.back(), temp, sumvar;
    vector<double> parameter_values(5), non_linear_parameters, ubar(n_models), ustarbar(n_models);
    NumericMatrix u(n_doses, n_models), ustar(n_doses, n_models);

    NumericMatrix opt_contrast(n_doses, n_models);

    sumvar = 0.0;
    for (j = 0; j < n_doses; j++) sumvar += 1.0 / var[j];

    // Set up dose-response models based on the initial values
    for (i = 0; i < n_models; i++) {

        non_linear_parameters = ExtractRow(user_specified, i);

        // Parameters of a standardized model
        parameter_values = StandardDRFunction(model_index[i], 0.0, 1.0, max_dose, non_linear_parameters);

        for (j = 0; j < n_doses; j++) {

          u(j, i) = DoseResponseFunction(dose_levels[j], model_index[i], parameter_values, 1.0);

        }

    }

    for (i = 0; i < n_models; i++) {

        ubar[i] = 0.0;
        for (j = 0; j < n_doses; j++) ubar[i] += u(j, i) / var[j];
        ubar[i] /= sumvar;    

        for (j = 0; j < n_doses; j++) ustar(j, i) = (u(j, i) - ubar[i]) / var[j];
        ustarbar[i] = 0;    
        for (j = 0; j < n_doses; j++) ustarbar[i] += ustar(j, i);
        ustarbar[i] /= (double) n_doses;
            
        temp = 0.0;
        for (j = 0; j < n_doses; j++) temp += Sq(ustar(j, i) - ustarbar[i]);  
        for (j = 0; j < n_doses; j++) opt_contrast(j, i) = (ustar(j, i) - ustarbar[i])/ sqrt(temp);            

    }                

    return opt_contrast;
}

//REM, 318
vector<double> FitLinearModel(const vector<double> &x, const vector<double> &y) {

    int i, n = x.size();
    vector<double> coef(2);
    double temp1 = 0.0, temp2 = 0.0;

    double x_mean = MeanVec(x);
    double y_mean = MeanVec(y);

    for(i = 0; i < n; i++) {
        temp1 += (x[i] - x_mean) * (y[i] - y_mean);
        temp2 += Sq(x[i] - x_mean);
    }

    coef[1] = temp1 / temp2;
    coef[0] = y_mean - coef[1] * x_mean;

    return coef;  
}

//REM, 338
class WLSFit: public MFuncGrad
{
private:
    const vector<double> X;
    const vector<double> Y;
    const NumericMatrix S;
    const vector<double> bounds;
    const int model;
public:
    WLSFit(const vector<double> x_, const vector<double> y_, const NumericMatrix s_, const vector<double> bounds_, const int model_) : X(x_), Y(y_), S(s_), bounds(bounds_), model(model_) {}

    double f_grad(Constvec& beta, Refvec grad)
    {

        int i, j, n = X.size(), n_parameters;
        double wls_criterion = 0.0, diffi, diffj, deni, denj, e0, e1, delta, ed50, emax;
        Eigen::VectorXd gradient; 

        if (model == 1) n_parameters = 2;
        if (model == 2 || model == 3) n_parameters = 3;
        if (model == 4) n_parameters = 4;

        gradient.resize(n_parameters); 
        for (i = 0; i < n_parameters; i++) gradient[i] = 0.0;           

        // Linear model
        if (model == 1) {    

            e0 = beta[0];    
            delta = beta[1];    

            for(i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    diffi = e0 + delta * X[i] - Y[i];                
                    diffj = e0 + delta * X[j] - Y[j];
                    wls_criterion += S(i, j) * diffi * diffj;
                    gradient[0] += S(i, j) * (diffi + diffj);
                    gradient[1] += S(i, j) * (X[j] * diffi + X[i] * diffj);
                }
            }
        }

        // Exponential model
        if (model == 2) {    

            e0 = beta[0];    
            e1 = beta[1];    
            delta = beta[2];
            if (delta < bounds[0]) delta = bounds[0];
            if (delta > bounds[1]) delta = bounds[1];

            for(i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    diffi = e0 + e1 * (exp(X[i] / delta) - 1.0) - Y[i];                
                    diffj = e0 + e1 * (exp(X[j] / delta) - 1.0) - Y[j];
                    wls_criterion += S(i, j) * diffi * diffj;
                    gradient[0] += S(i, j) * (diffi + diffj);
                    gradient[1] += S(i, j) * ((exp(X[i] / delta) - 1.0) * diffj + (exp(X[j] / delta) - 1.0) * diffi);
                    gradient[2] += S(i, j) * ((- e1 * exp(X[i] / delta) / Sq(delta)) * diffj + (- e1 * exp(X[j] / delta) / Sq(delta)) * diffi);
                }
            }

        }

        // Emax model
        if (model == 3) {    

            e0 = beta[0];    
            emax = beta[1];    
            ed50 = beta[2];    
            if (ed50 < bounds[0]) ed50 = bounds[0];
            if (ed50 > bounds[1]) ed50 = bounds[1];

            for(i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    diffi = e0 + emax * X[i] / (ed50 + X[i]) - Y[i];                
                    diffj = e0 + emax * X[j] / (ed50 + X[j]) - Y[j];
                    wls_criterion += S(i, j) * diffi * diffj;
                    gradient[0] += S(i, j) * (diffi + diffj);
                    gradient[1] += S(i, j) * (diffj * X[i] / (ed50 + X[i]) + diffi * X[j] / (ed50 + X[j]));
                    gradient[2] += S(i, j) * (diffj * (- emax * X[i] / Sq(ed50 + X[i])) + diffi * (- emax * X[j] / Sq(ed50 + X[j])));
                }
            }

        }

        // Logistic model
        if (model == 4) {    

            e0 = beta[0];    
            emax = beta[1];    
            ed50 = beta[2];    
            delta = beta[3];    
            if (ed50 < bounds[0]) ed50 = bounds[0]; 
            if (ed50 > bounds[1]) ed50 = bounds[1];
            if (delta < bounds[2]) delta = bounds[2];
            if (delta > bounds[3]) delta = bounds[3];

            for(i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    deni = 1.0 + exp((ed50 - X[i]) / delta);
                    denj = 1.0 + exp((ed50 - X[j]) / delta);
                    diffi = e0 + emax / deni - Y[i];                
                    diffj = e0 + emax / denj - Y[j];
                    wls_criterion += S(i, j) * diffi * diffj;
                    gradient[0] += S(i, j) * (diffi + diffj);
                    gradient[1] += S(i, j) * (diffj / deni + diffi / denj);
                    gradient[2] += S(i, j) * (diffj * (- emax * (deni - 1.0) / (delta * Sq(deni))) + diffi * (- emax * (denj - 1.0) / (delta * Sq(denj))));
                    gradient[3] += S(i, j) * (diffj * (emax * (deni - 1.0) * (ed50 - X[i]) / (Sq(delta * deni))) + diffi * (emax * (denj - 1.0) * (ed50 - X[j]) / (Sq(delta * denj))));
                }
            }

        }

        final_gradient = 0.0;
        for (i = 0; i < n_parameters; i++) final_gradient += abs(gradient[i]);

        const double f = wls_criterion;
        grad.noalias() = gradient;

        return f;
    }
};

//REM, 463
// Scaled inverse chi-square distribution
vector<double> ScaledInvChiSq(const int &n, const double &nu, const double &tau2) {
  
    int i;
    vector<double> result(n);

    for (i = 0; i < n; i++) 
        result[i] = 1.0 / Gamma(1, 
                                nu / 2.0, 
                                (nu * tau2) / 2.0)[0];

    return result;  

}

//REM, 477
vector<double> GeneratePosteriorSample(const int &nsamples, const double &mean_obs, const double &sd_obs, const int &n_obs, const double &nu_prior, const double &tau2_prior, const double &mean_prior, const double &k){
  
    int i; 

    double mean_post;
    vector<double> sample(nsamples), var_post(nsamples);

    // Variance follows and inverse chi-squared distribution
    var_post = ScaledInvChiSq(nsamples, 
                 nu_prior + n_obs, 
                 (nu_prior * tau2_prior + (n_obs - 1.0) * Sq(sd_obs) + (k * (n_obs + 0.0) * Sq(mean_prior - mean_obs)) / (k + n_obs)) / (nu_prior + n_obs + 0.0));

    for (i = 0; i < nsamples; i++) {
        mean_post = Normal(1, (k * mean_prior + (n_obs + 0.0) * mean_obs) / (k + n_obs + 0.0), sqrt(var_post[i] / (k + n_obs + 0.0)))[0];
        sample[i] = Normal(1, mean_post, sqrt(var_post[i]))[0];
    }

    return sample;

}

//REM, 498
void SetInitialValues(vector<ModelInformation> &model_information, const vector<double> &dose, const vector<double> &response, const double &max_dose, const vector<int> &model_index) {

    int i;
    vector<double> linear_coef(2);

    // Fit a linear regression model first    
    linear_coef = FitLinearModel(dose, response);

    double placebo_effect = linear_coef[0], max_effect = linear_coef[0] + linear_coef[1] * max_dose, d, d1, d2, p, p1, p2;

    // Initial values of the non-linear model parameters  
    vector<double> non_linear_parameters(3), temp_vec(1), bounds1(2), bounds2(4), temp_par;

    non_linear_parameters = fillvec(3, 0.0);

    temp_vec[0] = 0.0;

    ModelInformation current_model_information;

    for (i = 0; i < n_models; i++) {

        current_model_information.model_index = model_index[i];

        // Set initial values of the model parameters

        // Linear
        if (model_index[i] == 1) {

            current_model_information.n_parameters = 2;
            current_model_information.initial_values = linear_coef;
        }

        // Exponential
        if (model_index[i] == 2) {

            current_model_information.n_parameters = 3;

            non_linear_parameters[0] = max_dose;
            current_model_information.initial_values = ComputeDoseResponseFunctionParameters(2, placebo_effect, max_effect, max_dose, non_linear_parameters);

            // Use pre-defined bounds for non-linear model parameters
            bounds1[0] = 0.1 * max_dose;
            bounds1[1] = 2.0 * max_dose;

            current_model_information.bounds = bounds1;
        }

        // Emax
        if (model_index[i] == 3) {

            current_model_information.n_parameters = 3;

            d = 0.5 * max_dose;
            p = 0.5;
            non_linear_parameters[0] = d * (1.0 - p) / p;
            current_model_information.initial_values = ComputeDoseResponseFunctionParameters(3, placebo_effect, max_effect, max_dose, non_linear_parameters);

            // Use pre-defined bounds for non-linear model parameters
            bounds1[0] = 0.001 * max_dose;
            bounds1[1] = 1.5 * max_dose;

            current_model_information.bounds = bounds1;

        }

        // Logistic
        if (model_index[i] == 4) {

            current_model_information.n_parameters = 4;

            d1 = 0.33 * max_dose;
            p1 = 0.33;
            d2 = 0.66 * max_dose;
            p2 = 0.66;
            non_linear_parameters[0] = (d1 * Logit(p2) - d2 * Logit(p1)) / (Logit(p2) - Logit(p1));
            non_linear_parameters[1] = (d2 - d1)/(Logit(p2) - Logit(p1));
            current_model_information.initial_values = ComputeDoseResponseFunctionParameters(4, placebo_effect, max_effect, max_dose, non_linear_parameters);

            // Use pre-defined bounds for non-linear model parameters
            bounds2[0] = 0.001 * max_dose;
            bounds2[1] = 1.5 * max_dose;
            bounds2[2] = 0.01 * max_dose;
            bounds2[3] = 0.3 * max_dose;

            current_model_information.bounds = bounds2;

        }

        current_model_information.coef = fillvec(current_model_information.n_parameters, 0.0);

        current_model_information.status = -1;
        current_model_information.criterion = 10000.0;    
        current_model_information.target_dose = -1.0;
        current_model_information.convergence_criterion = -1.0;

        model_information[i] = current_model_information;        

    }

}

//REM, 599
// MCPMod analysis at the end of the trial
vector<double> MCPMod(const vector<int> &group_n, const vector<double> &group_mean, const vector<double> &group_sd, const vector<double> &dose_levels, const vector<int> &model_index, const vector<double> &non_linear_vector) {

    /*******************************************************************/

    double denom, numer, pooled_variance;

    int i, j, n_doses = group_n.size();

    vector<double> test_stat(n_models), sign_model(n_models), var(n_doses);

    NumericMatrix optimal_contrast, non_linear_matrix(4, 2);

    /*******************************************************************/

    // Non-linear parameters

    // Linear model
    non_linear_matrix(0, 0) = non_linear_vector[0]; 

    // Exponential model
    non_linear_matrix(1, 0) = non_linear_vector[1]; 

    // Emax model
    non_linear_matrix(2, 0) = non_linear_vector[2]; 

    // Logistic model
    non_linear_matrix(3, 0) = non_linear_vector[3]; 
    non_linear_matrix(3, 1) = non_linear_vector[4]; 

    for (j = 0; j < n_doses; j++) {
        var[j] = 1.0 / (group_n[j] + 0.0);
    }

    // Pooled variance
    pooled_variance = 0.0;
    for (j = 0; j < n_doses; j++) {
        pooled_variance += (group_n[j] - 1.0) * Sq(group_sd[j]);
    }
    pooled_variance /= (SumVecInt(group_n) - n_doses + 0.0);

    // Compute the optimal contrast 
    optimal_contrast = OptContrast(non_linear_matrix, model_index, direction_index_local, dose_levels, var);

    for (i = 0; i < n_models; i++) {
      denom = 0.0;
      numer = 0.0;
      for (j = 0; j < n_doses; j++) {
          numer = numer + optimal_contrast(j, i) * group_mean[j];
          denom = denom + Sq(optimal_contrast(j, i)) / (group_n[j] + 0.0);
      }
      test_stat[i] = numer / sqrt(pooled_variance * denom);
    }

    /*******************************************************************/

    return test_stat;
}

//REM, 659
// Fit dose-response models for the standard MCPMod method 
void FitDoseResponseModels(vector<ModelInformation> &model_information, const vector<double> &dose, const vector<double> &response, const NumericMatrix &cov_mat, const double &direction_index, const int &maxit, const double &max_gradient) {

    int i, j, n_parameters, convergence;

    vector<double> coef;

    double fopt, eps_f = 1e-08, eps_g = 1e-06;

    /*******************************************************************************/

    for (i = 0; i < n_models; i++) {

        n_parameters = model_information[i].n_parameters;
        Eigen::VectorXd beta(n_parameters);
        coef.clear();
        convergence = 1;

        for (j = 0; j < n_parameters; j++) {
            beta[j] = model_information[i].initial_values[j];               
        }

        WLSFit regression_fit(dose, response, cov_mat, model_information[i].bounds, model_information[i].model_index);
        model_information[i].status = optim_lbfgs(regression_fit, beta, fopt, maxit, eps_f, eps_g);

        // Criteria for determining convergence
        if (isnan(fopt) || isnan(final_gradient) || abs(final_gradient) > max_gradient || model_information[i].status < 0 || isnan(beta[0])) {
            convergence = 0;
            model_information[i].status = -1;
        }

        // Convergence
        if (convergence == 1) {

            // Apply constraints
            if (model_information[i].model_index == 2 || model_information[i].model_index == 3) {
                beta[2] = min(max(beta[2], model_information[i].bounds[0]), model_information[i].bounds[1]);           
            }     

            if (model_information[i].model_index == 4) {
                beta[2] = min(max(beta[2], model_information[i].bounds[0]), model_information[i].bounds[1]);                
                beta[3] = min(max(beta[3], model_information[i].bounds[2]), model_information[i].bounds[3]);                
            }

            for (j = 0; j < n_parameters; j++) coef.push_back(beta[j]);

            model_information[i].coef = coef;
            model_information[i].convergence_criterion = final_gradient;
            model_information[i].criterion = fopt + 2.0 * (n_parameters + 0.0); 

        }

    }
}


//REM, 714
vector<ModelInformation> ModelFit(const vector<int> &group_n, const vector<double> &group_mean, const vector<double> &group_sd, const vector<double> &dose_levels, const vector<int> &model_index, const vector<double> &non_linear_vector) {

    /*******************************************************************/

    double max_gradient = 1000.0, max_dose = dose_levels.back();

    int j, n_doses = group_n.size(), maxit = 50;

    NumericMatrix cov_mat(n_doses, n_doses), non_linear_matrix(4, 2);

    /*******************************************************************/

    // Non-linear parameters

    // Linear model
    non_linear_matrix(0, 0) = non_linear_vector[0]; 

    // Exponential model
    non_linear_matrix(1, 0) = non_linear_vector[1]; 

    // Emax model
    non_linear_matrix(2, 0) = non_linear_vector[2]; 

    // Logistic model
    non_linear_matrix(3, 0) = non_linear_vector[3]; 
    non_linear_matrix(3, 1) = non_linear_vector[4]; 

    for (j = 0; j < n_doses; j++) {
        cov_mat(j, j) = group_n[j];
    }

    // Vector of model information parameters
    vector<ModelInformation> model_information(n_models);

    SetInitialValues(model_information, dose_levels, group_mean, max_dose, model_index);

    FitDoseResponseModels(model_information, dose_levels, group_mean, cov_mat, direction_index_local, maxit, max_gradient);

    return model_information;

}

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
    int direction_index = as<int>(parameters["direction_index"]);

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

            if (endpoint_index == 1) test_result1 = TTest(control_sample, treatment_sample, 0.0, direction_index);
            if (endpoint_index == 2) test_result1 = PropTest(control_sample, treatment_sample, 0.0, direction_index);

            pvalue1 = test_result1.pvalue;

            trad_outcome = (pvalue1 <= alpha);        

            // Interim analysis 1      

            control_sample = ExtractSamples(overall_control_sample, 0, sample_size_ia1[0]);
            treatment_sample = ExtractSamples(overall_treatment_sample, 0, sample_size_ia1[1]);

            if (endpoint_index == 1) test_result1 = TTest(control_sample, treatment_sample, 0.0, direction_index);
            if (endpoint_index == 2) test_result1 = PropTest(control_sample, treatment_sample, 0.0, direction_index);

            // Conditional power
            cp[0] = CondPower(test_result1.test_stat, sample_size_ia1[0] + sample_size_ia1[1], sample_size_fa1[0] + sample_size_fa1[1], 1, 1, 0.0, alpha);

            outcome[0] = (cp[0] <= futility_threshold);

            // Interim analysis 2

            if (outcome[0] < 1.0) {

                control_sample = ExtractSamples(overall_control_sample, 0, sample_size_ia2[0]);
                treatment_sample = ExtractSamples(overall_treatment_sample, 0, sample_size_ia2[1]);

                if (endpoint_index == 1) test_result1 = TTest(control_sample, treatment_sample, 0.0, direction_index);
                if (endpoint_index == 2) test_result1 = PropTest(control_sample, treatment_sample, 0.0, direction_index);

                // Conditional power
                cp[1] = CondPower(test_result1.test_stat, sample_size_ia2[0] + sample_size_ia2[1], sample_size_fa1[0] + sample_size_fa1[1], 1, 1, 0.0, alpha);

                pvalue1 = test_result1.pvalue;

            }

            if (outcome[0] < 1.0 && outcome[1] < 1.0) outcome[2] = (cp[1] > underpowered_zone[0] && cp[1] <= underpowered_zone[1]); 

            // Final analysis (before sample size adjustment)
            if (outcome[0] < 1.0 && outcome[1] < 1.0 && outcome[2] < 1.0) {

                control_sample = ExtractSamples(overall_control_sample, 0, sample_size_fa1[0]);
                treatment_sample = ExtractSamples(overall_treatment_sample, 0, sample_size_fa1[1]);

                if (endpoint_index == 1) test_result2 = TTest(control_sample, treatment_sample, 0.0, direction_index);
                if (endpoint_index == 2) test_result2 = PropTest(control_sample, treatment_sample, 0.0, direction_index);

                pvalue2 = test_result2.pvalue;

                outcome[3] = (pvalue2 <= alpha);        

            } 

            // Final analysis (after sample size adjustment)
            if (outcome[0] < 1.0 && outcome[1] < 1.0 && outcome[2] > 0.0) {

                // Updated sample size in both arms
                sample_size_increase = UpdatedEventCount(test_result1.test_stat, sample_size_ia2[0] + sample_size_ia2[1], sample_size_fa1[0] + sample_size_fa1[1], sample_size_fa2[0] + sample_size_fa2[1], 1, 1, 0.0, target_power, alpha); 

                control_sample = ExtractSamples(overall_control_sample, sample_size_ia2[0], sample_size_ia2[0] + (int) (ratio * sample_size_increase));
                treatment_sample = ExtractSamples(overall_treatment_sample, sample_size_ia2[1], sample_size_ia2[1] + (int) ((1.0 - ratio) * sample_size_increase));

                if (endpoint_index == 1) test_result2 = TTest(control_sample, treatment_sample, 0.0, direction_index);
                if (endpoint_index == 2) test_result2 = PropTest(control_sample, treatment_sample, 0.0, direction_index);

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

            test_result1 = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

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

            test_result1 = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);
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

                test_result1 = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);
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

                test_result2 = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);
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

                test_result2 = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);
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

    int treatment_count = as<int>(parameters["treatment_count"]);
    int mult_test_index = as<int>(parameters["mult_test_index"]);
    vector<double> weight = as<vector<double>>(parameters["weight"]);
    vector<double> transition = as<vector<double>>(parameters["transition"]);

    vector<double> dropout_parameter = as<vector<double>>(parameters["dropout_parameter"]);

    int endpoint_index = as<int>(parameters["endpoint_index"]);
    int direction_index = as<int>(parameters["direction_index"]);

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

    double cp_threshold, ad_outcome;

    int i, j, sim, current_stratum, ntotal, narms = sample_size.size();

    vector<double> adj_pvalue, select_flag(narms - 1), futility_flag(narms - 1);

    ntotal = SumVecInt(sample_size);

    AllSurvivalData survival_data;

    OutcomeCensor outcome_censor, outcome_censor_control, outcome_censor_treatment;

    TestResult test_result;

    vector<double> look_time(3), event_count(3), dropout, hr(narms), cp1(narms - 1), cp2(narms - 1), sortedcp2(narms - 1), pvalue(narms - 1), overall_control_sample, control_sample, treatment_sample, temp_vec, pvalue_trad(narms - 1), trad_outcome(narms - 1);

    vector<int> stratum_list, sample_list;

    NumericMatrix sim_results(nsims, 3 * narms - 2), overall_treatment_sample(narms - 1, max_sample_size); 

    /*******************************************************************/

    for (sim = 0; sim < nsims; sim++) {

        hr = fillvec(narms, 0.0);
        cp1 = fillvec(narms - 1, 0.0);
        cp2 = fillvec(narms - 1, 0.0);
        pvalue = fillvec(narms - 1, 1.0);
        trad_outcome = fillvec(narms - 1, 0.0);
        select_flag = fillvec(narms - 1, 0.0);
        futility_flag = fillvec(narms - 1, 0.0);

        look_time = fillvec(3, 0.0);
        event_count = fillvec(3, 0.0);
        pvalue_trad = fillvec(narms - 1, 0.0);

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

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);

                // Conditional power
                cp1[i] = CondPower(test_result.test_stat, sample_size_ia1[0] + sample_size_ia1[i + 1], sample_size_fa[0] + sample_size_fa[i + 1], 1, 1, 0.0, alpha);

                // Futility stopping rule
                futility_flag[i] = (cp1[i] <= futility_threshold);

            }

            // Interim analysis 2

            control_sample = ExtractSamples(overall_control_sample, 0, sample_size_ia2[0]);

            for (i = 0; i < narms - 1; i++) {

                treatment_sample.clear();
                for (j = 0; j < sample_size_ia2[i + 1]; j++) treatment_sample.push_back(overall_treatment_sample(i, j));

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);

                // Conditional power
                cp2[i] = CondPower(test_result.test_stat, sample_size_ia2[0] + sample_size_ia2[i + 1], sample_size_fa[0] + sample_size_fa[i + 1], 1, 1, 0.0, alpha);

            }

            for (i = 0; i < narms - 1; i++) sortedcp2[i] = cp2[i];

            // Treatment selection rule
            sort(sortedcp2.begin(), sortedcp2.end());

            cp_threshold = sortedcp2[narms - treatment_count - 1];

            for (i = 0; i < narms - 1; i++) {

                if (futility_flag[i] < 1.0 && cp2[i] >= cp_threshold) {
                    select_flag[i] = 1.0;    
                }

            }

            // Final analysis

            for (i = 0; i < narms - 1; i++) {

                control_sample = ExtractSamples(overall_control_sample, 0, sample_size_fa[0]);

                if (select_flag[i] > 0.0) {

                    treatment_sample.clear();
                    for (j = 0; j < sample_size_fa[i + 1]; j++) treatment_sample.push_back(overall_treatment_sample(i, j));

                    if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
                    if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);

                    pvalue[i] = test_result.pvalue;
                
                } 

            }

            // Apply the multiplicity adjustment

            adj_pvalue = TradMultAdj(mult_test_index, pvalue, weight, transition);

           for (i = 0; i < narms - 1; i++) {

              if (adj_pvalue[i] < alpha) ad_outcome = 1.0;         

           }

            // Traditional analysis

            control_sample = ExtractSamples(overall_control_sample, 0, sample_size_fa[0]);

            for (i = 0; i < narms - 1; i++) {

                treatment_sample.clear();
                for (j = 0; j < sample_size_fa[i + 1]; j++) treatment_sample.push_back(overall_treatment_sample(i, j));

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);

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
                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

                // Conditional power
                cp1[i] = CondPower(test_result.test_stat, event_count_ia1, event_count_fa, 1, 1, 0.0, alpha);

                // Futility stopping rule
                futility_flag[i] = (cp1[i] <= futility_threshold);

            }

            // Interim analysis 2

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
                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

                // Conditional power
                cp2[i] = CondPower(test_result.test_stat, event_count_ia2, event_count_fa, 1, 1, 0.0, alpha);

            }

            // Treatment selection rule

            for (i = 0; i < narms - 1; i++) sortedcp2[i] = cp2[i];

            // Treatment selection rule
            sort(sortedcp2.begin(), sortedcp2.end());

            cp_threshold = sortedcp2[narms - treatment_count - 1];

            for (i = 0; i < narms - 1; i++) {

                if (futility_flag[i] < 1.0 && cp2[i] >= cp_threshold) {
                    select_flag[i] = 1.0;    
                }

            }

            // Final analysis

            stratum_list.clear();
            for (i = 0; i < narms; i++) stratum_list.push_back(i);        
            look_time[2] = FindMilestone(stratum_list, survival_data.stratum, survival_data.os_local_start, event_count_fa);

            for (i = 0; i < narms - 1; i++) {

                // Control
                stratum_list.clear();
                stratum_list.push_back(0);
                outcome_censor_control = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);

                if (select_flag[i] > 0.0) {

                    // Treatment
                    stratum_list.clear();
                    stratum_list.push_back(i + 1);
                    outcome_censor_treatment = ExtractOutcomeCensor(stratum_list, survival_data.stratum, survival_data.start, survival_data.os_local, survival_data.os_local_censor, look_time[2]);
        
                    test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

                    pvalue[i] = test_result.pvalue;

                }
            }

            // Apply the multiplicity adjustment

            adj_pvalue = TradMultAdj(mult_test_index, pvalue, weight, transition);

           for (i = 0; i < narms - 1; i++) {

              if (adj_pvalue[i] < alpha) ad_outcome = 1.0;         

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
    
                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

                pvalue_trad[i] = test_result.pvalue;

                trad_outcome[i] = (pvalue_trad[i] <= alpha && futility_flag[i] < 1.0);        

            }

        }
        // Time-to-event endpoints

        // Summary statistics for each simulation run

        sim_results(sim, 0) = ad_outcome;
        for (i = 0; i < narms - 1; i++) sim_results(sim, 1 + i) = trad_outcome[i];
        for (i = 0; i < narms - 1; i++) sim_results(sim, 1 + narms - 1 + i) = futility_flag[i];
        for (i = 0; i < narms - 1; i++) sim_results(sim, 1 + 2 * narms - 2 + i) = select_flag[i];

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
    int direction_index = as<int>(parameters["direction_index"]);

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

    NumericMatrix sim_results(nsims, 15); 

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

            if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
            if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);

            // Conditional power
            cp = CondPower(test_result.test_stat, SumVecInt(sample_size_ia1), SumVecInt(sample_size_fa), 1, 1, 0.0, alpha);

            // Futility stopping rule
            futility_flag = (cp <= futility_threshold);

            // Interim analysis 2

            // Biomarker-negative subpopulation

            control_sample = ExtractSamples(sample0, 0, sample_size_ia2[0]);
            treatment_sample = ExtractSamples(sample2, 0, sample_size_ia2[2]);

            effect_size_minus = ComputeEffectSize(control_sample, treatment_sample, endpoint_index, 1, direction_index);

            // Biomarker-positive subpopulation

            control_sample = ExtractSamples(sample1, 0, sample_size_ia2[1]);
            treatment_sample = ExtractSamples(sample3, 0, sample_size_ia2[3]);

            effect_size_plus = ComputeEffectSize(control_sample, treatment_sample, endpoint_index, 1, direction_index);

            hypothesis_selection_outcome = HypothesisSelection(effect_size_minus, effect_size_plus, influence, interaction); 

            // Final analysis

            // All patients are included in the final analysis
            if (hypothesis_selection_outcome[0] > 0.0 || hypothesis_selection_outcome[2] > 0.0) {

                // Overall population

                control_sample = CombineVec(ExtractSamples(sample0, 0, sample_size_fa[0]), ExtractSamples(sample1, 0, sample_size_fa[1]));
                treatment_sample = CombineVec(ExtractSamples(sample2, 0, sample_size_fa[2]), ExtractSamples(sample3, 0, sample_size_fa[3]));

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);

                ad_outcome[0] = (test_result.pvalue <= alpha && futility_flag < 1.0);

                // Biomarker-positive subpopulation

                control_sample = ExtractSamples(sample1, 0, sample_size_fa[1]);
                treatment_sample = ExtractSamples(sample3, 0, sample_size_fa[3]);

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);

                ad_outcome[1] = (test_result.pvalue <= alpha && futility_flag < 1.0);

            } 

            // Only biomarker-positive patients are included in the final analysis
            if (hypothesis_selection_outcome[1] > 0.0) {

                control_sample = ExtractSamples(sample1, 0, sample_size_fa[1]);
                treatment_sample = ExtractSamples(sample3, 0, sample_size_fa[3]);

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

                ad_outcome[1] = (test_result.pvalue <= alpha && futility_flag < 1.0);

            } 

            // Traditional analysis

            // Overall population

            control_sample = CombineVec(ExtractSamples(sample0, 0, sample_size_fa[0]), ExtractSamples(sample1, 0, sample_size_fa[1]));
            treatment_sample = CombineVec(ExtractSamples(sample2, 0, sample_size_fa[2]), ExtractSamples(sample3, 0, sample_size_fa[3]));

            if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
            if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);

            trad_outcome[0] = (test_result.pvalue <= alpha && futility_flag < 1.0);        
            // Biomarker-positive subpopulation

            control_sample = ExtractSamples(sample1, 0, sample_size_fa[1]);
            treatment_sample = ExtractSamples(sample3, 0, sample_size_fa[3]);

            if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
            if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);

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
            test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);
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

            if (direction_index == 1) {
                effect_size_minus = -log(hr_minus);
                effect_size_plus = -log(hr_plus);
            }

            if (direction_index == 2) {
                effect_size_minus = log(hr_minus);
                effect_size_plus = log(hr_plus);
            }

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

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

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

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

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

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

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

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

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

            test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

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

            test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

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
    int direction_index = as<int>(parameters["direction_index"]);

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
                    test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
                }

                // Binary endpoint
                if (endpoint_index == 2) {
                    treatment_sample = Binary(interim[i + 1], rates[i + 1]);
                    test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);
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

                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

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

// [[Rcpp::export]]
List ADRandC(const List &parameters_arg) {

    // TODO: Get this variables from parameters?
    //int n_models = 4, direction_index_local = 1;

    List parameters(parameters_arg);

    vector<int> stage_sample_size = as<vector<int>>(parameters["stage_sample_size"]);
    vector<double> mean = as<vector<double>>(parameters["mean"]);
    vector<double> sd = as<vector<double>>(parameters["sd"]);
    vector<double> dose_levels = as<vector<double>>(parameters["dose_levels"]);

    vector<int> model_index = as<vector<int>>(parameters["model_index"]);
    vector<double> non_linear_vector = as<vector<double>>(parameters["non_linear_vector"]);

    int direction_index = as<int>(parameters["direction_index"]);

    double treatment_period = as<double>(parameters["treatment_period"]);    

    double enrollment_period = as<double>(parameters["enrollment_period"]);    
    double enrollment_parameter = as<double>(parameters["enrollment_parameter"]);    
    double dropout_rate = as<double>(parameters["dropout_rate"]);    

    double delta = as<double>(parameters["delta"]);    
    double ratio_placebo = as<double>(parameters["ratio_placebo"]);    
    double balance = as<double>(parameters["balance"]);    

    int nsims = as<int>(parameters["nsims"]);    
    int n_per_arm = as<int>(parameters["n_per_arm"]);    

    /*******************************************************************/

    double denom, denominator, current_criterion, planned_end, overrun1, overrun2;

    int i, j, sim, stage, n_doses = dose_levels.size(), n_total = SumVecInt(stage_sample_size), n_stages = stage_sample_size.size(), nsamples = 1000, current_sample_size;

    vector<int> n(n_doses), sample_n(n_doses);
    vector<double> traditional(n_models), adaptive(n_models), sample_mean(n_doses), predicted_mean(n_doses), sample_sd(n_doses), patient_start, randomization_ratio(n_doses), post_prob(n_doses - 1), sample_sum(n_doses), sample_sum_sq(n_doses), placebo_sample, current_sample, criterion(n_models), linear_coef(2), exp_coef(3), emax_coef(3), logistic_coef(4), model_weight(n_models), coef;

    NumericMatrix sample_n_matrix(nsims, n_doses), traditional_matrix(nsims, n_models),adaptive_matrix(nsims, n_models), prediction(n_models, n_doses), stage_sample_size_matrix(nsims, n_stages);

    vector<ModelInformation> model_information;
    
    /*******************************************************************/

    for (sim = 0; sim < nsims; sim++) { 

        // Traditional design

        sample_n = FillVecInt(n_doses, 0);
        sample_sum = fillvec(n_doses, 0.0);
        sample_sum_sq = fillvec(n_doses, 0.0);

        // Generate the data and compute descriptive statistics
        for (i = 0; i < n_doses; i++) {
            n[i] = (int) round(n_per_arm * (1.0 - dropout_rate));
            if (direction_index == 1) current_sample = Normal(n[i], mean[i], sd[i]); 
            if (direction_index == 2) current_sample = Normal(n[i], -mean[i], sd[i]);
            sample_sum[i] = sum(current_sample);
            sample_sum_sq[i] = sumsq(current_sample);
            sample_mean[i] = sample_sum[i] / (double) n[i];
            sample_sd[i] = sqrt(sample_sum_sq[i] / (double) n[i] - Sq(sample_mean[i]));
        }        

        // Perform the MCPMod analysis at the final analysis in the traditional design
        traditional = MCPMod(n, sample_mean, sample_sd, dose_levels, model_index, non_linear_vector);

        // Adaptive design

        sample_n = FillVecInt(n_doses, 0);
        sample_sum = fillvec(n_doses, 0.0);
        sample_sum_sq = fillvec(n_doses, 0.0);

        // Patient enrollment times
        patient_start = Enrollment(n_total, enrollment_period, 2, enrollment_parameter); 
        sort(patient_start.begin(), patient_start.end());

        // Compute the initial randomization ratios
        randomization_ratio[0] = ratio_placebo;
        for (i = 1; i < n_doses; i++) {
            randomization_ratio[i] = (1.0 - ratio_placebo) / (n_doses - 1.0);
        }

        overrun1 = 0.0;
        current_sample_size = 0;

        // Actual number of patients in each stage 
        for (stage = 0; stage < n_stages; stage++) {

            current_sample_size += stage_sample_size[stage];

            // Planned end of current stage
            planned_end = patient_start[current_sample_size - 1];

            // Number of overrun patients enrolled between the planned end and interim look
            overrun2 = 0.0;
            if (stage < n_stages - 1) {

                if (patient_start[current_sample_size] < planned_end + treatment_period) {
                    j = 0;
                    while(current_sample_size + j < n_total && patient_start[current_sample_size + j] < planned_end + treatment_period) {
                        j++;
                    }
                    overrun2 = (double) j;
                }
            } 

            stage_sample_size_matrix(sim, stage) = stage_sample_size[stage] + overrun2 - overrun1;
            overrun1 = overrun2;

        }

        for (stage = 0; stage < n_stages; stage++) {

            // Generate the data and compute descriptive statistics
            for (i = 0; i < n_doses; i++) {
                n[i] = (int) round(stage_sample_size_matrix(sim, stage) * randomization_ratio[i] * (1.0 - dropout_rate));
                if (n[i] <= 0) n[i] = 1;
                if (direction_index == 1) current_sample = Normal(n[i], mean[i], sd[i]); 
                if (direction_index == 2) current_sample = Normal(n[i], -mean[i], sd[i]);
                sample_n[i] += n[i];
                sample_sum[i] += sum(current_sample);
                sample_sum_sq[i] += sumsq(current_sample);
                sample_mean[i] = sample_sum[i] / (double) sample_n[i];
                sample_sd[i] = sqrt(sample_sum_sq[i] / (double) sample_n[i] - Sq(sample_mean[i]));
            }

            // Fit the candidate dose-response model, apply model averaging and compute predicted means

            model_information = ModelFit(sample_n, sample_mean, sample_sd, dose_levels, model_index, non_linear_vector);

            for (i = 0; i < 2; i++) linear_coef[i] = model_information[0].coef[i];
            for (i = 0; i < 3; i++) exp_coef[i] = model_information[1].coef[i];
            for (i = 0; i < 3; i++) emax_coef[i] = model_information[2].coef[i];
            for (i = 0; i < 4; i++) logistic_coef[i] = model_information[3].coef[i];

            // Model-specific predicted means
            for (i = 0; i < n_models; i++) {    

                if (i == 0) coef = linear_coef;
                if (i == 1) coef = exp_coef;
                if (i == 2) coef = emax_coef;
                if (i == 3) coef = logistic_coef;

                for (j = 0; j < n_doses; j++) prediction(i, j) = DoseResponseFunction(dose_levels[j], model_index[i], coef, 1); 

            }

            // Model weights
            for (i = 0; i < n_models; i++) {

                model_weight[i] = 0.0;

                if (model_information[i].status >= 0) {

                    current_criterion = model_information[i].criterion;
                    denominator = 0.0;

                    for (j = 0; j < n_models; j++) {

                        if (model_information[j].status >= 0) denominator += exp(- 0.5 * (model_information[j].criterion - current_criterion));

                    }

                    if (abs(denominator) > 0.0001) model_weight[i] = 1.0 / denominator;

                }

            }

            for (j = 0; j < n_doses; j++) {

                predicted_mean[j] = 0.0;
                for (i = 0; i < n_models; i++) predicted_mean[j] += (model_weight[i] * prediction(i, j));

            }

            // Compute the posterior probabilities of target efficacy

            // Sample from the posterior distribution to create a placebo sample
            placebo_sample = GeneratePosteriorSample(nsamples, predicted_mean[0], sample_sd[0], sample_n[0], -1.0, 0.0, 0.0, 0.0);

            for (i = 1; i < n_doses; i++) {
                current_sample = GeneratePosteriorSample(nsamples, predicted_mean[i], sample_sd[i], sample_n[i], -1.0, 0.0, 0.0, 0.0);
                post_prob[i - 1] = 0.0;
                for (j = 0; j < nsamples; j++) post_prob[i - 1] += (current_sample[j] >= placebo_sample[j] + delta);
                post_prob[i - 1] /= (nsamples + 0.0);   

                if (isnan(post_prob[i - 1])) post_prob[i - 1] = 1.0;
            }

            denom = 0.0;
            for (i = 0; i < n_doses - 1; i++) denom += pow(post_prob[i], balance);

            // Update the randomization ratio
            randomization_ratio[0] = ratio_placebo;

            if (denom > 0.0) {
                for (i = 1; i < n_doses; i++) {
                    randomization_ratio[i] = (1.0 - ratio_placebo) * pow(post_prob[i - 1], balance) / denom;
                }
            } else {
                for (i = 1; i < n_doses; i++) {
                    randomization_ratio[i] = (1.0 - ratio_placebo) / (n_doses - 1.0);
                }               
            }

        }
        // End of stages

        // Perform the MCPMod analysis at the final analysis in the adaptive design
        adaptive = MCPMod(sample_n, sample_mean, sample_sd, dose_levels, model_index, non_linear_vector);

        // Save the results
        for (i = 0; i < n_doses; i++) {
            sample_n_matrix(sim, i) = sample_n[i];
        }

        for (i = 0; i < n_models; i++) traditional_matrix(sim, i) = traditional[i];
        for (i = 0; i < n_models; i++) adaptive_matrix(sim, i) = adaptive[i];

    }

    /*******************************************************************/

    // Return simulation results

    return List::create(Named("n") = sample_n_matrix,
                        Named("traditional") = traditional_matrix,
                        Named("adaptive") = adaptive_matrix,
                        Named("stage_sample_size") = stage_sample_size_matrix);

}
// [[Rcpp::export]]
List MultAdjC1(const List &parameters_arg) {

    List parameters(parameters_arg);

    vector<int> sample_size = as<vector<int>>(parameters["sample_size_adj"]);
    int max_sample_size = as<int>(parameters["max_sample_size"]);

    int ncomparisons = as<int>(parameters["n_comparisons"]);
    int mult_test_index = as<int>(parameters["mult_test_index"]);

    int direction = as<int>(parameters["direction_index"]);

    int endpoint_index = as<int>(parameters["endpoint_index"]);

    vector<double> means = as<vector<double>>(parameters["means"]);
    vector<double> sds = as<vector<double>>(parameters["sds"]);
    vector<double> rates = as<vector<double>>(parameters["rates"]);

    vector<double> weight = as<vector<double>>(parameters["weights"]);
    vector<double> transition = as<vector<double>>(parameters["ctransition"]);

    int nsims = as<int>(parameters["nsims"]);

    /*******************************************************************/

    vector<double> temp_vec, control_sample, treatment_sample, pvalue(ncomparisons), adj_pvalue(ncomparisons);

    vector<int> stratum_list;

    int i, j, sim;

    NumericMatrix sim_results(nsims, ncomparisons + ncomparisons), overall_treatment_sample(ncomparisons, max_sample_size); 

    /*******************************************************************/

    for (sim = 0; sim < nsims; sim++) {

        // Normal and binary endpoints
        if (endpoint_index == 1 || endpoint_index == 2) {

            // Control samples
            if (endpoint_index == 1) control_sample = Normal(sample_size[0], means[0], sds[0]);
            if (endpoint_index == 2) control_sample = Binary(sample_size[0], rates[0]);

            // Treatment samples
            for (i = 0; i < ncomparisons; i++) {

                for (j = 0; j < max_sample_size; j++) overall_treatment_sample(i, j) = 0.0;

                if (endpoint_index == 1) temp_vec = Normal(sample_size[i + 1], means[i + 1], sds[i + 1]);
                if (endpoint_index == 2) temp_vec = Binary(sample_size[i + 1], rates[i + 1]);

                for (j = 0; j < sample_size[i + 1]; j++) overall_treatment_sample(i, j) = temp_vec[j];

            }

            // Final analysis

            for (i = 0; i < ncomparisons; i++) {

                treatment_sample.clear();
                for (j = 0; j < sample_size[i + 1]; j++) treatment_sample.push_back(overall_treatment_sample(i, j));

                if (endpoint_index == 1) pvalue[i] = TTest(control_sample, treatment_sample, 0.0, direction).pvalue;
                if (endpoint_index == 2) pvalue[i] = PropTest(control_sample, treatment_sample, 0.0, direction).pvalue;

            }

        }
        // Normal and binary endpoints

        // Multiplicity adjustment
        adj_pvalue = TradMultAdj(mult_test_index, pvalue, weight, transition);

        // Summary statistics for each simulation run

        for (i = 0; i < ncomparisons; i++) sim_results(sim, i) = pvalue[i];
        for (i = 0; i < ncomparisons; i++) sim_results(sim, i + ncomparisons) = adj_pvalue[i];

    }
    // End of simulations

    /*******************************************************************/

    // Return simulation results

    return List::create(Named("sim_results") = sim_results);

}
// End of MultAdjC1

// [[Rcpp::export]]
List MultAdjC2(const List &parameters_arg) {

    List parameters(parameters_arg);

    vector<int> sample_size = as<vector<int>>(parameters["sample_size_adj"]);

    int nendpoints = as<int>(parameters["n_endpoints"]);
    int mult_test_index = as<int>(parameters["mult_test_index"]);

    int direction = as<int>(parameters["direction_index"]);

    int endpoint_index = as<int>(parameters["endpoint_index"]);

    NumericMatrix corr = as<NumericMatrix>(parameters["endpoint_correlation"]);
    double corr_sum = as<double>(parameters["corr_sum"]);

    vector<double> control_mean = as<vector<double>>(parameters["control_mean"]);
    vector<double> treatment_mean = as<vector<double>>(parameters["treatment_mean"]);
    vector<double> control_sd = as<vector<double>>(parameters["control_sd"]);
    vector<double> treatment_sd = as<vector<double>>(parameters["treatment_sd"]);
    vector<double> control_rate = as<vector<double>>(parameters["control_rate"]);
    vector<double> treatment_rate = as<vector<double>>(parameters["treatment_rate"]);

    vector<double> weight = as<vector<double>>(parameters["weights"]);
    vector<double> transition = as<vector<double>>(parameters["ctransition"]);

    int nsims = as<int>(parameters["nsims"]);

    /*******************************************************************/

    double overall_pvalue, overall_test_stat;

    vector<double> temp_vec, control_sample, treatment_sample, pvalue(nendpoints), test_stat(nendpoints), adj_pvalue(nendpoints);

    vector<int> stratum_list;

    int i, j, sim, ntotal;

    ntotal = SumVecInt(sample_size);

    TestResult test_result;

    NumericMatrix sim_results(nsims, nendpoints + nendpoints), overall_control_sample(nendpoints, sample_size[0]), overall_treatment_sample(nendpoints, sample_size[1]); 

    /*******************************************************************/

    for (sim = 0; sim < nsims; sim++) {

        // Normal endpoints
        if (endpoint_index == 1) {

            // Control sample
            for (i = 0; i < sample_size[0]; i++) {
                temp_vec = MVNormal(nendpoints, control_mean, control_sd, corr); 
                for (j = 0; j < nendpoints; j++) overall_control_sample(j, i) = temp_vec[j];  
            }      

            // Treatment sample
            for (i = 0; i < sample_size[1]; i++) {
                temp_vec = MVNormal(nendpoints, treatment_mean, treatment_sd, corr); 
                for (j = 0; j < nendpoints; j++) overall_treatment_sample(j, i) = temp_vec[j];  
            }      

            // Final analysis
            for (i = 0; i < nendpoints; i++) {

                control_sample.clear();
                for (j = 0; j < sample_size[0]; j++) control_sample.push_back(overall_control_sample(i, j));

                treatment_sample.clear();
                for (j = 0; j < sample_size[1]; j++) treatment_sample.push_back(overall_treatment_sample(i, j));

                test_result = TTest(control_sample, treatment_sample, 0.0, direction);    
                pvalue[i] = test_result.pvalue;
                test_stat[i] = test_result.test_stat;

            }

        }
        // Normal endpoints

        // Binary endpoints
        if (endpoint_index == 2) {

            // Control sample
            for (i = 0; i < sample_size[0]; i++) {
                temp_vec = MVNormal(nendpoints, control_mean, control_sd, corr); 
                for (j = 0; j < nendpoints; j++) {
                    // Convert to binary outcomes
                    overall_control_sample(j, i) = (temp_vec[j] <= rcpp_qnorm(control_rate[j]));  
                }
            }      

            // Treatment sample
            for (i = 0; i < sample_size[1]; i++) {
                temp_vec = MVNormal(nendpoints, treatment_mean, treatment_sd, corr); 
                for (j = 0; j < nendpoints; j++) {
                    // Convert to binary outcomes
                    overall_treatment_sample(j, i) = (temp_vec[j] <= rcpp_qnorm(treatment_rate[j]));  
                }
            }      

            // Final analysis
            for (i = 0; i < nendpoints; i++) {

                control_sample.clear();
                for (j = 0; j < sample_size[0]; j++) control_sample.push_back(overall_control_sample(i, j));

                treatment_sample.clear();
                for (j = 0; j < sample_size[1]; j++) treatment_sample.push_back(overall_treatment_sample(i, j));

                test_result = PropTest(control_sample, treatment_sample, 0.0, direction);    
                pvalue[i] = test_result.pvalue;
                test_stat[i] = test_result.test_stat;

            }

        }
        // Binary endpoints

        // Multiplicity adjustment
        if (mult_test_index <= 6) {

            adj_pvalue = TradMultAdj(mult_test_index, pvalue, weight, transition);

            // Summary statistics for each simulation run

            for (i = 0; i < nendpoints; i++) sim_results(sim, i) = pvalue[i];
            for (i = 0; i < nendpoints; i++) sim_results(sim, i + nendpoints) = adj_pvalue[i];

        }

        // Global test
        if (mult_test_index == 7) {

            overall_test_stat = sum(test_stat) / sqrt(corr_sum);
            overall_pvalue = 1.0 - rcpp_pt(overall_test_stat, 0.5 * (ntotal - 2.0) * (1.0 + 1.0 / Sq(nendpoints + 0.0)));

            // Summary statistics for each simulation run
            for (i = 0; i < nendpoints; i++) sim_results(sim, i) = pvalue[i];
            sim_results(sim, nendpoints) = overall_pvalue;
            for (i = 1; i < nendpoints; i++) sim_results(sim, i + nendpoints) = 0.0;

        }

    }
    // End of simulations

    /*******************************************************************/

    // Return simulation results

    return List::create(Named("sim_results") = sim_results);

}
// End of MultAdjC2

// [[Rcpp::export]]
List MultAdjC3(const List &parameters_arg) {

    List parameters(parameters_arg);

    vector<int> sample_size = as<vector<int>>(parameters["sample_size_adj"]);
    int max_sample_size = as<int>(parameters["max_sample_size"]);

    int ncomparisons = as<int>(parameters["n_comparisons"]);
    int nendpoints = as<int>(parameters["n_endpoints"]);
    int mult_test_index = as<int>(parameters["mult_test_index"]);

    int direction = as<int>(parameters["direction_index"]);

    int endpoint_index = as<int>(parameters["endpoint_index"]);

    NumericMatrix corr = as<NumericMatrix>(parameters["endpoint_correlation"]);

    vector<double> control_mean = as<vector<double>>(parameters["control_mean"]);
    vector<double> control_sd = as<vector<double>>(parameters["control_sd"]);
    NumericMatrix treatment_mean = as<NumericMatrix>(parameters["treatment_mean"]);
    NumericMatrix treatment_sd = as<NumericMatrix>(parameters["treatment_sd"]);
    vector<double> control_rate = as<vector<double>>(parameters["control_rate"]);
    NumericMatrix treatment_rate = as<NumericMatrix>(parameters["treatment_rate"]);

    int nsims = as<int>(parameters["nsims"]);

    vector<double> gamma = as<vector<double>>(parameters["mult_test_gamma"]);
    int method = as<int>(parameters["mult_method_index"]);

    /*******************************************************************/

    int i, j, k, l, sim, nhypotheses = nendpoints * ncomparisons;

    vector<double> temp_vec(nendpoints), control_sample, treatment_sample, pvalue(nhypotheses), adj_pvalue(nhypotheses), current_treatment_mean, current_treatment_sd, temp_mean, temp_sd;

    vector<int> stratum_list;

    NumericMatrix sim_results(nsims, nhypotheses + nhypotheses), overall_control_sample(nendpoints, sample_size[0]), overall_treatment_sample(nendpoints, max_sample_size); 

    VariableCuboid overall_treatment_sample_mat(ncomparisons, nendpoints, max_sample_size);

    /*******************************************************************/

    for (sim = 0; sim < nsims; sim++) {

        // Normal endpoints
        if (endpoint_index == 1) {

            // Control sample
            for (i = 0; i < sample_size[0]; i++) {
                temp_vec = MVNormal(nendpoints, control_mean, control_sd, corr); 
                for (j = 0; j < nendpoints; j++) overall_control_sample(j, i) = temp_vec[j];  
            }      

            // Treatment samples
            for (k = 0; k < ncomparisons; k++) {

                current_treatment_mean = ExtractColumn(treatment_mean, k);
                current_treatment_sd = ExtractColumn(treatment_sd, k);

                for (i = 0; i < sample_size[k + 1]; i++) {

                    temp_vec = MVNormal(nendpoints, current_treatment_mean, current_treatment_sd, corr); 

                    for (j = 0; j < nendpoints; j++) {
                        overall_treatment_sample_mat(k, j, i) = temp_vec[j];  
                    }

                }      

            }

            // Hypothesis tests
            l = 0;
            for (i = 0; i < nendpoints; i++) {

                control_sample.clear();
                for (j = 0; j < sample_size[0]; j++) control_sample.push_back(overall_control_sample(i, j));

                for (k = 0; k < ncomparisons; k++) {

                    treatment_sample.clear();
                    for (j = 0; j < sample_size[k + 1]; j++) treatment_sample.push_back(overall_treatment_sample_mat(k, i, j)); 

                    pvalue[l] = TTest(control_sample, treatment_sample, 0.0, direction).pvalue;

                    l++;

                }

            }

        }
        // Normal endpoints
        
        // Binary endpoints
        if (endpoint_index == 2) {

            temp_mean = fillvec(nendpoints, 0.0);
            temp_sd = fillvec(nendpoints, 1.0);                

            // Control sample
            for (i = 0; i < sample_size[0]; i++) {

                temp_vec = MVNormal(nendpoints, temp_mean, temp_sd, corr); 
                // Convert to binary outcomes
                for (j = 0; j < nendpoints; j++) overall_control_sample(j, i) = (temp_vec[j] <= rcpp_qnorm(control_rate[j]));  
            }      

            // Treatment samples
            for (k = 0; k < ncomparisons; k++) {

                for (i = 0; i < sample_size[k + 1]; i++) {

                    temp_vec = MVNormal(nendpoints, temp_mean, temp_sd, corr); 
                    // Convert to binary outcomes
                    for (j = 0; j < nendpoints; j++) overall_treatment_sample_mat(k, j, i) = (temp_vec[j] <= rcpp_qnorm(treatment_rate(j, k)));   

                }      

            }

            // Hypothesis tests
            l = 0;
            for (i = 0; i < nendpoints; i++) {

                control_sample.clear();
                for (j = 0; j < sample_size[0]; j++) control_sample.push_back(overall_control_sample(i, j));

                for (k = 0; k < ncomparisons; k++) {

                    treatment_sample.clear();
                    for (j = 0; j < sample_size[k + 1]; j++) treatment_sample.push_back(overall_treatment_sample_mat(k, i, j)); 

                    pvalue[l] = PropTest(control_sample, treatment_sample, 0.0, direction).pvalue;

                    l++;

                }

            }

        }
        // Binary endpoints

        // Multiplicity adjustment
        adj_pvalue = MixtureProcAdjP(nendpoints, ncomparisons, pvalue, mult_test_index, gamma, method);

        // Summary statistics for each simulation run

        for (i = 0; i < nhypotheses; i++) sim_results(sim, i) = pvalue[i];
        for (i = 0; i < nhypotheses; i++) sim_results(sim, i + nhypotheses) = adj_pvalue[i];


    }
    // End of simulations

    /*******************************************************************/

    // Return simulation results

    return List::create(Named("sim_results") = sim_results);

}
// End of MultAdjC3

// Two vectors of simulation results
struct SimulationResults {
    vector<double> pval;
    vector<double> eff;
};

// Convert a vector
NumericVector FromVectorXd(const VectorXd &eigen_vec) {

    int i, a = eigen_vec.size();  

    NumericVector res(a);

    for (i = 0; i < a; i++) res[i] = eigen_vec[i];    

    return(res);

}

// Convert a matrix
NumericMatrix FromMatrixXd(const MatrixXd &eigen_mat) {

    int i, j, a = eigen_mat.rows();  

    NumericMatrix res(a, a);

    for (i = 0; i < a; i++) {
        for (j = 0; j < a; j++) {
            res(i, j) = eigen_mat(i, j);    
        }
    }

    return(res);

}

// Convert a matrix
MatrixXd ToMatrixXd(const NumericMatrix &mat) {

    int i, j, a = mat.nrow();  

    MatrixXd res(a, a);

    for (i = 0; i < a; i++) {
        for (j = 0; j < a; j++) {
            res(i, j) = mat(i, j);    
        }
    }

    return(res);

}


// # nocov start
// Extract a subset of a matrix
NumericMatrix ExtractMat(const NumericMatrix &mat, const vector<int> &id, const int &value) {

    int i, j, k, c = 0, a = mat.nrow(), b = mat.ncol(), n = id.size();
    for (i = 0; i < n; i++) {
        if (id[i] == value) c++;
    }
    
    NumericMatrix res(c, b);

    k = 0;
    for (i = 0; i < a; i++) {
        if (id[i] == value) {
            for (j = 0; j < b; j++) {
                res(k, j) = mat(i, j);
                k++;
            }
        }
    }

    return res;

}
// # nocov end

// Transpose a matrix
NumericMatrix TransMat(const NumericMatrix &mat) {

    int i, j, a = mat.nrow(), b = mat.ncol();
    
    NumericMatrix res(b, a);

    for (i = 0; i < a; i++) {
       for (j = 0; j < b; j++) {
            res(j, i) = mat(i, j); 
       }
    }   

    return res;

}    

// Multiply two matrices
NumericMatrix MultMat(const NumericMatrix &mat1, const NumericMatrix &mat2) {

    int i, j, k, a = mat1.nrow(), b = mat1.ncol(), c = mat2.ncol();  
    double sum;
    NumericMatrix res(a, c);

    for (i = 0; i < a; i++) {
        for (j = 0; j < c; j++) {
            sum = 0.0;
            for (k = 0; k < b; k++) sum += mat1(i, k) * mat2(k, j);
            res(i, j) = sum;    
        }
    }

    return res;

}

// Add two matrices
NumericMatrix AddMat(const NumericMatrix &mat1, const NumericMatrix &mat2) {

    int i, j, a = mat1.nrow(), b = mat1.ncol();  
    NumericMatrix res(a, b);

    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++) {
            res(i, j) = mat1(i, j) + mat2(i, j);    
        }
    }

    return res;

}

// Subtract two matrices
NumericMatrix SubtractMat(const NumericMatrix &mat1, const NumericMatrix &mat2) {

    int i, j, a = mat1.nrow(), b = mat1.ncol();  
    NumericMatrix res(a, b);

    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++) {
            res(i, j) = mat1(i, j) - mat2(i, j);    
        }
    }

    return res;

}

// Invert a square matrix
NumericMatrix InvMat(const NumericMatrix &mat) {

    int a = mat.nrow();  

    MatrixXd eigen_mat(a, a), inv_eigen_mat(a, a);
    NumericMatrix res(a, a);

    eigen_mat = ToMatrixXd(mat); 

    inv_eigen_mat = eigen_mat.inverse();

    res = FromMatrixXd(inv_eigen_mat);

    return res;

}

// Square root of a square matrix
NumericMatrix InvSqRootMat(const NumericMatrix &mat) {

    int i, j, a = mat.nrow();
    
    NumericMatrix res(a, a);
    MatrixXd eigen_mat(a, a);

    eigen_mat = ToMatrixXd(mat);

    // Singular value decomposition
    JacobiSVD<MatrixXd>svd(eigen_mat, Eigen::ComputeThinU|Eigen::ComputeThinV);

    NumericMatrix u, v, sqrt_diag(a, a), temp1;
    NumericVector s;

    u = FromMatrixXd(svd.matrixU());
    v = FromMatrixXd(svd.matrixV());
    s = FromVectorXd(svd.singularValues());

    // Square roots of the diagonal elements
    for (i = 0; i < a; i++) {
        for (j = 0; j < a; j++) {
            if (i == j) sqrt_diag(i, j) = 1.0 / sqrt(s(i)); else sqrt_diag(i, j) = 0.0;   
        }
    }    

    // Square root of the original matrix
    temp1 = MultMat(u, sqrt_diag);
    res = MultMat(temp1, TransMat(v));

    return res;

}

// Generate a vector of random cluster sizes
vector<int> RandomClusterSize(const int &sample_size, const vector<double> &proportion) {

    int i, j, m = proportion.size();
    vector<int> cluster_size(m);
    vector<double> unif;

    unif = Uniform(sample_size, 0.0, 1.0);

    for (i = 0; i < sample_size; i++) {

        if (unif[i] <= proportion[0]) {

            cluster_size[0]++;

        } else {

            if (unif[i] > proportion[m - 2]) {

                cluster_size[m - 1]++;

            } else {

                for (j = 0; j < m - 2; j++) {

                    if (unif[i] > proportion[j] && unif[i] <= proportion[j + 1]) {
                        cluster_size[j + 1]++;
                    }
                }
            }
        }
    }

    return(cluster_size);

}

// Compute initial values of the binary GEE model
NumericMatrix BinInitValues(const vector<double> &y, const NumericMatrix &x) {

    int i, j, n = y.size();
    vector<double> p(n), var(n);
    NumericMatrix xvar(n, 2), r(n, 1), outcome(n, 1), est(2, 1), inc(2, 1), pr, temp1, temp2;

    for (i = 0; i < n; i++) outcome(i, 0) = 2.0 * y[i] - 1.0;

    temp1 = InvMat(MultMat(TransMat(x), x));
    temp2 = MultMat(temp1, TransMat(x));
    est = MultMat(temp2, outcome);

    for (j = 0; j < 2; j++) {

        pr = MultMat(x, est);

        for (i = 0; i < n; i++) {
            p[i] = 1.0 / (1.0 + exp(-pr(i, 0)));
            var[i] = p[i] * (1.0 - p[i]);
            r(i, 0) = y[i] - p[i];
            xvar(i, 0) = x(i, 0) * var[i];
            xvar(i, 1) = x(i, 1) * var[i];
        }

        outcome = MultMat(TransMat(x), r);
        temp1 = InvMat(MultMat(TransMat(x), xvar));
        inc = MultMat(temp1, outcome);

        est = AddMat(est, inc);

    }

    return est;

}

SimulationResults ContGEE(const vector<double> &y, const NumericMatrix &x, const vector<int> &id, const vector<int> &id_freq, const int &direction_index, const int &max_iterations, const double &min_increment) {

    int i, j, k, cycle, pat, nobs = y.size(), npats = id_freq.size(), npars = 2, npoints, continue_flag, convergence_flag;

    vector<double> est(npars), stderr_sandwich(npars), stderr_kc(npars), stderr_md(npars), pval(4), eff(3);

    double corr, ave, increment, corr_off_diagonal, corr_diagonal, pvalue_sandwich, pvalue_kc, pvalue_md;

    NumericMatrix next_coef(npars, 1), current_coef(npars, 1), est_mat(npars, 1), hessian_mat(npars, npars), gradient_mat_sandwich(npars, npars), gradient_mat_kc(npars, npars), gradient_mat_md(npars, npars), pred, temp1, temp2, temp3, temp4, temp5, inv_corr_mat, cov_mat_sandwich, cov_mat_kc, cov_mat_md, inv_hessian_mat;

    SimulationResults sim_results;

    // Initialize the values    
    ave = 0.0;    
    for (i = 0; i < npats; i++) ave += 0.5 * (id_freq[i] + 0.0) * (id_freq[i] - 1.0);
    corr = 0;
    next_coef = FillMat(next_coef, 0.0);

    // Stopping condition for model fitting
    continue_flag = 1;
    cycle = 0;

    convergence_flag = 1;

    while (continue_flag == 1) {

        // Coefficient estimates in the current step
        for (i = 0; i < npars; i++) current_coef(i, 0) = next_coef(i, 0);

        est_mat = FillMat(est_mat, 0.0);
        hessian_mat = FillMat(hessian_mat, 0.0);
        gradient_mat_sandwich = FillMat(gradient_mat_sandwich, 0.0);
        gradient_mat_kc = FillMat(gradient_mat_kc, 0.0);
        gradient_mat_md = FillMat(gradient_mat_md, 0.0);

        corr_off_diagonal = 0.0;
        corr_diagonal = 0.0;

        // Loop over patients/clusters
        for (pat = 0; pat < npats; pat++) {

            npoints = id_freq[pat];

            // Correlation matrix for the current patient (under the assumption of an exchangeable working correlation matrix)
            NumericMatrix corr_mat(npoints, npoints);
            for (i = 0; i < npoints; i++) {
                for (j = 0; j < npoints; j++) {
                    if (i == j) corr_mat(i, j) = 1.0; else corr_mat(i, j) = corr;    
                }
            }            
            inv_corr_mat = InvMat(corr_mat);

            // Extract the x and y values for the current patient/cluster 
            vector<double> y_cluster(npoints);
            NumericMatrix x_cluster(npoints, npars);
            k = 0;
            for (i = 0; i < nobs; i++) {
                if (id[i] == pat + 1) {
                    y_cluster[k] = y[i];
                    for (j = 0; j < npars; j++) x_cluster(k, j) = x(i, j);                   
                    k++;    
                }
            } 

            // Predicted means
            NumericMatrix residual(npoints, 1);
            pred = MultMat(x_cluster, current_coef);

            for (i = 0; i < npoints; i++) residual(i, 0) = y_cluster[i] - pred(i, 0);
            for (i = 0; i < npoints; i++) corr_diagonal += Sq(residual(i, 0));    
            for (j = 0; j < npoints - 1; j++) {
                for (k = j + 1; k < npoints; k++) corr_off_diagonal += residual(j, 0) * residual(k, 0);   
            }
 
            // Update the estimate matrix
            temp1 = MultMat(inv_corr_mat, residual);
            temp2 = MultMat(TransMat(x_cluster), temp1);
            est_mat = AddMat(est_mat, temp2);

            // Update the Hessian matrix
            temp1 = MultMat(inv_corr_mat, x_cluster);
            temp2 = MultMat(TransMat(x_cluster), temp1);
            hessian_mat = AddMat(hessian_mat, temp2);

            // Update the gradient matrix (uncorrected sandwich estimator)
            temp1 = MultMat(TransMat(x_cluster), inv_corr_mat);
            temp2 = MultMat(temp1, residual);
            temp3 = MultMat(temp2, TransMat(temp2));
            gradient_mat_sandwich = AddMat(gradient_mat_sandwich, temp3);

        }

        // Simple model-based covariance estimate
        inv_hessian_mat = InvMat(hessian_mat);

        // Separate loop over patients/clusters to compute bias-corrected covariance estimators
        for (pat = 0; pat < npats; pat++) {

            npoints = id_freq[pat];

            NumericMatrix identity_mat(npoints, npoints), leverage_mat(npoints, npoints);

            // Correlation matrix for the current patient (under the assumption of an exchangeable working correlation matrix)
            NumericMatrix corr_mat(npoints, npoints);
            for (i = 0; i < npoints; i++) {
                for (j = 0; j < npoints; j++) {
                    if (i == j) {
                        corr_mat(i, j) = 1.0;
                        identity_mat(i, j) = 1.0;
                    }
                    else {
                        corr_mat(i, j) = corr;   
                        identity_mat(i, j) = 0.0;
                    } 
                }
            }            
            inv_corr_mat = InvMat(corr_mat);

            // Extract the x and y values for the current patient/cluster 
            vector<double> y_cluster(npoints);
            NumericMatrix x_cluster(npoints, npars);
            k = 0;
            for (i = 0; i < nobs; i++) {
                if (id[i] == pat + 1) {
                    y_cluster[k] = y[i];
                    for (j = 0; j < npars; j++) x_cluster(k, j) = x(i, j);                   
                    k++;    
                }
            } 

            temp1 = MultMat(x_cluster, inv_hessian_mat);
            temp2 = MultMat(temp1, TransMat(x_cluster));
            leverage_mat = MultMat(temp2, inv_corr_mat);

            // Predicted means
            NumericMatrix residual(npoints, 1);
            pred = MultMat(x_cluster, current_coef);
            for (i = 0; i < npoints; i++) residual(i, 0) = y_cluster[i] - pred(i, 0);

            // Update the gradient matrix (sandwich estimator by Kauermann and Carroll)
            temp1 = MultMat(TransMat(x_cluster), inv_corr_mat);
            temp2 = InvSqRootMat(SubtractMat(identity_mat, leverage_mat)); 

            temp3 = MultMat(temp1, temp2);
            temp4 = MultMat(temp3, residual);
            temp5 = MultMat(temp4, TransMat(temp4));
            gradient_mat_kc = AddMat(gradient_mat_kc, temp5);

            // Update the gradient matrix (sandwich estimator by Mancl and DeRouen)
            temp1 = MultMat(TransMat(x_cluster), inv_corr_mat);
            temp2 = InvMat(SubtractMat(identity_mat, leverage_mat));

            temp3 = MultMat(temp1, temp2);
            temp4 = MultMat(temp3, residual);
            temp5 = MultMat(temp4, TransMat(temp4));
            gradient_mat_md = AddMat(gradient_mat_md, temp5);


        }

        // Update the common correlation coefficient (under the assumption of an exchangeable working correlation matrix)
        corr = corr_off_diagonal * (npats - npars + 0.0) / (corr_diagonal * (ave - npars + 0.0));

        // Update the coefficient estimates
        temp1 = MultMat(inv_hessian_mat, est_mat);
        next_coef = AddMat(current_coef, temp1);

        // Covariance matrix for the coefficient estimates (uncorrected sandwich estimator)
        temp1 = MultMat(gradient_mat_sandwich, inv_hessian_mat);
        cov_mat_sandwich = MultMat(inv_hessian_mat, temp1);

        // Covariance matrix for the coefficient estimates (sandwich estimator by Kauermann and Carroll)
        temp1 = MultMat(gradient_mat_kc, inv_hessian_mat);
        cov_mat_kc = MultMat(inv_hessian_mat, temp1);

        // Covariance matrix for the coefficient estimates (sandwich estimator by Mancl and DeRouen)
        temp1 = MultMat(gradient_mat_md, inv_hessian_mat);
        cov_mat_md = MultMat(inv_hessian_mat, temp1);

        cycle++;

        increment = 0;
        for (i = 0; i < npars; i++) increment += abs(next_coef[i] - current_coef[i]);

        // Stopping condition for model fitting
        if (cycle >= max_iterations || increment <= min_increment) continue_flag = 0;

    }
    // End of model fitting       

    if (cycle >= max_iterations) convergence_flag = 0;
  
    // Model parameter estimates and standard errors
    for (i = 0; i < npars; i++) {
        est[i] = next_coef(i, 0);
        stderr_sandwich[i] = sqrt(cov_mat_sandwich(i, i));
        stderr_kc[i] = sqrt(cov_mat_kc(i, i));
        stderr_md[i] = sqrt(cov_mat_md(i, i));
    }

    // One-sided p-values for the treatment effect
    pvalue_sandwich = 1.0 - rcpp_pnorm(est[1] / stderr_sandwich[1]);
    pvalue_kc = 1.0 - rcpp_pnorm(est[1] / stderr_kc[1]);
    pvalue_md = 1.0 - rcpp_pnorm(est[1] / stderr_md[1]);

    pval[0] = convergence_flag;
    pval[1] = pvalue_sandwich;
    pval[2] = pvalue_kc;
    pval[3] = pvalue_md;

    eff[0] = convergence_flag;
    eff[1] = est[0];
    eff[2] = est[1];

    sim_results.pval = pval;
    sim_results.eff = eff;

    return sim_results;

}
// End of ContGEE

SimulationResults BinGEE(const vector<double> &y, const NumericMatrix &x, const vector<int> &id, const vector<int> &id_freq, const int &direction_index, const int &max_iterations, const double &min_increment) {

    int i, j, k, cycle, pat, nobs = y.size(), npats = id_freq.size(), npars = 2, npoints, continue_flag, convergence_flag;

    vector<double> est(npars), stderr_sandwich(npars), stderr_kc(npars), stderr_md(npars), pval(4), eff(3);

    double corr, ave, increment, corr_off_diagonal, corr_diagonal, pvalue_sandwich, pvalue_kc, pvalue_md;

    NumericMatrix next_coef(npars, 1), current_coef(npars, 1), est_mat(npars, 1), hessian_mat(npars, npars), gradient_mat_sandwich(npars, npars), gradient_mat_kc(npars, npars), gradient_mat_md(npars, npars), pred, temp1, temp2, temp3, temp4, temp5, inv_corr_mat, cov_mat_sandwich, cov_mat_kc, cov_mat_md, inv_hessian_mat;

    SimulationResults sim_results;

    // Initialize the values    
    ave = 0.0;    
    for (i = 0; i < npats; i++) ave += 0.5 * (id_freq[i] + 0.0) * (id_freq[i] - 1.0);
    corr = 0;
    // Compute initial values of the binary GEE model
    next_coef = BinInitValues(y, x);

    // Stopping condition for model fitting
    continue_flag = 1;
    cycle = 0;

    convergence_flag = 1;

    while (continue_flag == 1) {

        // Coefficient estimates in the current step
        for (i = 0; i < npars; i++) current_coef(i, 0) = next_coef(i, 0);

        est_mat = FillMat(est_mat, 0.0);
        hessian_mat = FillMat(hessian_mat, 0.0);
        gradient_mat_sandwich = FillMat(gradient_mat_sandwich, 0.0);
        gradient_mat_kc = FillMat(gradient_mat_kc, 0.0);
        gradient_mat_md = FillMat(gradient_mat_md, 0.0);

        corr_off_diagonal = 0.0;
        corr_diagonal = 0.0;

        // Loop over patients/clusters
        for (pat = 0; pat < npats; pat++) {

            npoints = id_freq[pat];

            // Extract the x and y values for the current patient/cluster 
            vector<double> y_cluster(npoints);
            NumericMatrix x_cluster(npoints, npars);
            k = 0;
            for (i = 0; i < nobs; i++) {
                if (id[i] == pat + 1) {
                    y_cluster[k] = y[i];
                    for (j = 0; j < npars; j++) x_cluster(k, j) = x(i, j);                   
                    k++;    
                }
            } 

            // Predicted proportions
            pred = MultMat(x_cluster, current_coef);
            vector<double> prop(npoints);
            for (i = 0; i < npoints; i++) prop[i] = 1.0 / (1.0 + exp(-pred(i, 0)));

            // Correlation matrix for the current patient (under the assumption of an exchangeable working correlation matrix)
            NumericMatrix corr_mat(npoints, npoints);
            for (i = 0; i < npoints; i++) {
                for (j = 0; j < npoints; j++) {
                    if (i == j) corr_mat(i, j) = 1.0; else corr_mat(i, j) = corr;    
                }
            }            
 
            // Variance matrices
            NumericMatrix sqrt_var_mat(npoints, npoints);
            for (i = 0; i < npoints; i++) {
                for (j = 0; j < npoints; j++) {
                    if (i == j) sqrt_var_mat(i, j) = sqrt(prop[i] * (1.0 - prop[i])); else sqrt_var_mat(i, j) = 0.0;    
                } 
            } 
            temp1 = MultMat(sqrt_var_mat, corr_mat);
            corr_mat = MultMat(temp1, sqrt_var_mat);
            inv_corr_mat = InvMat(corr_mat);

            // Design matrix
            NumericMatrix design_mat(npoints, npars);
            for (i = 0; i < npoints; i++) {
                design_mat(i, 0) = prop[i] * (1.0 - prop[i]) * x_cluster(i, 0);
                design_mat(i, 1) = prop[i] * (1.0 - prop[i]) * x_cluster(i, 1);
            } 

            NumericMatrix residual(npoints, 1);
            for (i = 0; i < npoints; i++) residual(i, 0) = y_cluster[i] - prop[i];
            for (i = 0; i < npoints; i++) corr_diagonal += Sq(residual(i, 0));    
            for (j = 0; j < npoints - 1; j++) {
                for (k = j + 1; k < npoints; k++) corr_off_diagonal += residual(j, 0) * residual(k, 0);   
            }
 
            // Update the estimate matrix
            temp1 = MultMat(inv_corr_mat, residual);
            temp2 = MultMat(TransMat(design_mat), temp1);
            est_mat = AddMat(est_mat, temp2);

            // Update the Hessian matrix
            temp1 = MultMat(inv_corr_mat, design_mat);
            temp2 = MultMat(TransMat(design_mat), temp1);
            hessian_mat = AddMat(hessian_mat, temp2);

            // Update the gradient matrix (uncorrected sandwich estimator)
            temp1 = MultMat(TransMat(design_mat), inv_corr_mat);
            temp2 = MultMat(temp1, residual);
            temp3 = MultMat(temp2, TransMat(temp2));
            gradient_mat_sandwich = AddMat(gradient_mat_sandwich, temp3);

        }

        // Simple model-based covariance estimate
        inv_hessian_mat = InvMat(hessian_mat);

        // Separate loop over patients/clusters to compute bias-corrected covariance estimators
        for (pat = 0; pat < npats; pat++) {

            npoints = id_freq[pat];

            NumericMatrix identity_mat(npoints, npoints), leverage_mat(npoints, npoints);

            // Extract the x and y values for the current patient/cluster 
            vector<double> y_cluster(npoints);
            NumericMatrix x_cluster(npoints, npars);
            k = 0;
            for (i = 0; i < nobs; i++) {
                if (id[i] == pat + 1) {
                    y_cluster[k] = y[i];
                    for (j = 0; j < npars; j++) x_cluster(k, j) = x(i, j);                   
                    k++;    
                }
            } 

            // Predicted proportions
            pred = MultMat(x_cluster, current_coef);
            vector<double> prop(npoints);
            for (i = 0; i < npoints; i++) prop[i] = 1.0 / (1.0 + exp(-pred(i, 0)));

            // Correlation matrix for the current patient (under the assumption of an exchangeable working correlation matrix)
            NumericMatrix corr_mat(npoints, npoints);
            for (i = 0; i < npoints; i++) {
                for (j = 0; j < npoints; j++) {
                    if (i == j) {
                        corr_mat(i, j) = 1.0;
                        identity_mat(i, j) = 1.0;
                    }
                    else {
                        corr_mat(i, j) = corr;   
                        identity_mat(i, j) = 0.0;
                    } 
                }
            }            

            // Variance matrices
            NumericMatrix sqrt_var_mat(npoints, npoints);
            for (i = 0; i < npoints; i++) {
                for (j = 0; j < npoints; j++) {
                    if (i == j) sqrt_var_mat(i, j) = sqrt(prop[i] * (1.0 - prop[i])); else sqrt_var_mat(i, j) = 0.0;    
                } 
            } 
            temp1 = MultMat(sqrt_var_mat, corr_mat);
            corr_mat = MultMat(temp1, sqrt_var_mat);
            inv_corr_mat = InvMat(corr_mat);

            // Design matrix
            NumericMatrix design_mat(npoints, npars);
            for (i = 0; i < npoints; i++) {
                design_mat(i, 0) = prop[i] * (1.0 - prop[i]) * x_cluster(i, 0);
                design_mat(i, 1) = prop[i] * (1.0 - prop[i]) * x_cluster(i, 1);
            } 

            NumericMatrix residual(npoints, 1);
            for (i = 0; i < npoints; i++) residual(i, 0) = y_cluster[i] - prop[i];

            temp1 = MultMat(design_mat, inv_hessian_mat);
            temp2 = MultMat(temp1, TransMat(design_mat));
            leverage_mat = MultMat(temp2, inv_corr_mat);

            // Update the gradient matrix (sandwich estimator by Kauermann and Carroll)
            temp1 = MultMat(TransMat(design_mat), inv_corr_mat);
            temp2 = InvSqRootMat(SubtractMat(identity_mat, leverage_mat)); 
            temp3 = MultMat(temp1, temp2);
            temp4 = MultMat(temp3, residual);
            temp5 = MultMat(temp4, TransMat(temp4));
            gradient_mat_kc = AddMat(gradient_mat_kc, temp5);

            // Update the gradient matrix (sandwich estimator by Mancl and DeRouen)
            temp1 = MultMat(TransMat(design_mat), inv_corr_mat);
            temp2 = InvMat(SubtractMat(identity_mat, leverage_mat));
            temp3 = MultMat(temp1, temp2);
            temp4 = MultMat(temp3, residual);
            temp5 = MultMat(temp4, TransMat(temp4));
            gradient_mat_md = AddMat(gradient_mat_md, temp5);

        }

        // Update the common correlation coefficient (under the assumption of an exchangeable working correlation matrix)
        corr = corr_off_diagonal * (npats - npars + 0.0) / (corr_diagonal * (ave - npars + 0.0));

        // Update the coefficient estimates
        temp1 = MultMat(inv_hessian_mat, est_mat);
        next_coef = AddMat(current_coef, temp1);

        // Covariance matrix for the coefficient estimates (uncorrected sandwich estimator)
        temp1 = MultMat(gradient_mat_sandwich, inv_hessian_mat);
        cov_mat_sandwich = MultMat(inv_hessian_mat, temp1);

        // Covariance matrix for the coefficient estimates (sandwich estimator by Kauermann and Carroll)
        temp1 = MultMat(gradient_mat_kc, inv_hessian_mat);
        cov_mat_kc = MultMat(inv_hessian_mat, temp1);

        // Covariance matrix for the coefficient estimates (sandwich estimator by Mancl and DeRouen)
        temp1 = MultMat(gradient_mat_md, inv_hessian_mat);
        cov_mat_md = MultMat(inv_hessian_mat, temp1);

        cycle++;

        increment = 0;
        for (i = 0; i < npars; i++) increment += abs(next_coef[i] - current_coef[i]);

        // Stopping condition for model fitting
        if (cycle >= max_iterations || increment <= min_increment) continue_flag = 0;

    }
    // End of model fitting       

    if (cycle >= max_iterations) convergence_flag = 0;

    // Model parameter estimates and standard errors
    for (i = 0; i < npars; i++) {
        est[i] = next_coef(i, 0);
        stderr_sandwich[i] = sqrt(cov_mat_sandwich(i, i));
        stderr_kc[i] = sqrt(cov_mat_kc(i, i));
        stderr_md[i] = sqrt(cov_mat_md(i, i));
    }

    // One-sided p-values for the treatment effect
    pvalue_sandwich = 1.0 - rcpp_pnorm(est[1] / stderr_sandwich[1]);
    pvalue_kc = 1.0 - rcpp_pnorm(est[1] / stderr_kc[1]);
    pvalue_md = 1.0 - rcpp_pnorm(est[1] / stderr_md[1]);

    /*******************************************************************/

    pval[0] = convergence_flag;
    pval[1] = pvalue_sandwich;
    pval[2] = pvalue_kc;
    pval[3] = pvalue_md;

    eff[0] = convergence_flag;
    eff[1] = est[0];
    eff[2] = est[1];

    sim_results.pval = pval;
    sim_results.eff = eff;

    return sim_results;

}
// End of BinGEE

// [[Rcpp::export]]
List ClustRandGEEC(const List &parameters_arg) {

    List parameters(parameters_arg);

    int narms = as<int>(parameters["narms"]);
    int endpoint_index = as<int>(parameters["endpoint_index"]);
    int cluster_index = as<int>(parameters["cluster_index"]);
    vector<int> control_cluster_size_arg = as<vector<int>>(parameters["control_cluster_size"]);
    vector<double> control_cluster_cum = as<vector<double>>(parameters["control_cluster_cum"]);
    vector<int> treatment_cluster_size_vector = as<vector<int>>(parameters["treatment_cluster_size_vector"]);
    IntegerMatrix treatment_cluster_size_matrix = as<IntegerMatrix>(parameters["treatment_cluster_size_matrix"]);
    vector<double> treatment_cluster_cum_vector = as<vector<double>>(parameters["treatment_cluster_cum_vector"]);
    NumericMatrix treatment_cluster_cum_matrix = as<NumericMatrix>(parameters["treatment_cluster_cum_matrix"]);

    vector<int> sample_size = as<vector<int>>(parameters["sample_size"]);

    vector<double> means = as<vector<double>>(parameters["means"]);
    vector<double> within_cluster_sds = as<vector<double>>(parameters["within_cluster_sds"]);
    vector<double> between_cluster_sds = as<vector<double>>(parameters["between_cluster_sds"]);

    int direction_index = as<int>(parameters["direction_index"]);

    int mult_test_index = as<int>(parameters["mult_test_index"]);
    vector<double> weight = as<vector<double>>(parameters["weight"]);
    vector<double> transition = as<vector<double>>(parameters["transition"]);

    double control_alpha = as<double>(parameters["control_alpha"]);
    double control_beta = as<double>(parameters["control_beta"]);
    vector<double> treatment_alpha = as<vector<double>>(parameters["treatment_alpha"]);
    vector<double> treatment_beta = as<vector<double>>(parameters["treatment_beta"]);

    int max_iterations = as<int>(parameters["max_iterations"]);
    double min_increment = as<double>(parameters["min_increment"]);

    int nsims = as<int>(parameters["nsims_per_core"]);

    /*******************************************************************/

    double cluster_term;

    int i, j, k, sim, control_size, treatment_size;

    NumericMatrix pval_results(nsims, 4 * (narms - 1)), coef_results(nsims, 3 * (narms - 1)), pval(narms - 1, 4), eff(narms - 1, 3); 

    vector<int> id, id_control, id_freq, id_freq_control, temp_vec_int, control_cluster_size, treatment_cluster_size;

    vector<double> y, y_control, x, x_control, control_sample, current_treatment_sample, temp_vec, pvalue(narms - 1), adj_pvalue, cluster_size_results;

    SimulationResults results;

    /*******************************************************************/

    cluster_size_results.clear();

    for (sim = 0; sim < nsims; sim++) {

        // Fixed cluster sizes
        if (cluster_index == 1) {

            control_cluster_size.clear();
            for (i = 0; i < control_cluster_size_arg.size(); i++) control_cluster_size.push_back(control_cluster_size_arg[i]);

        }

        // Random cluster sizes
        if (cluster_index == 2) {

            // Generate a vector of random cluster sizes
            control_cluster_size = RandomClusterSize(sample_size[0], control_cluster_cum);

            cluster_size_results.insert(cluster_size_results.end(), control_cluster_size.begin(), control_cluster_size.end());

        }

        control_size = control_cluster_size.size();

        // Generate data for the common control arm
        x_control.clear();
        y_control.clear();
        id_control.clear();

        // Number of patients in each cluster
        id_freq_control.clear();
        id_freq_control.insert(id_freq_control.end(), control_cluster_size.begin(), control_cluster_size.end());

        for (i = 0; i < control_size; i++) {

            // Normal endpoint
            if (endpoint_index == 1) {

                cluster_term = Normal(1, 0.0, between_cluster_sds[0])[0];
                control_sample = Normal(control_cluster_size[i], means[0] + cluster_term, within_cluster_sds[0]);
            }

            // Binary endpoint
            if (endpoint_index == 2) {

                cluster_term = Beta(1, control_alpha, control_beta)[0];
                control_sample = Binary(control_cluster_size[i], cluster_term);

            }

            temp_vec = fillvec(control_cluster_size[i], 0.0);
            x_control.insert(x_control.end(), temp_vec.begin(), temp_vec.end());

            temp_vec = control_sample;
            y_control.insert(y_control.end(), temp_vec.begin(), temp_vec.end());

            temp_vec_int = FillVecInt(control_cluster_size[i], i + 1);
            id_control.insert(id_control.end(), temp_vec_int.begin(), temp_vec_int.end());

        }

        // Generate data for each treatment arm
        for (k = 1; k < narms; k++) {

            // Fixed cluster sizes
            if (cluster_index == 1) {

                treatment_cluster_size.clear();
                if (narms == 2) {
                
                    for (i = 0; i < treatment_cluster_size_vector.size(); i++) treatment_cluster_size.push_back(treatment_cluster_size_vector[i]); 

                } else {

                    for (i = 0; i < treatment_cluster_size_matrix.ncol(); i++) treatment_cluster_size.push_back(treatment_cluster_size_matrix(k - 1, i));   // # nocov

                }

            }

            // Random cluster sizes
            if (cluster_index == 2) {

                if (narms == 2) {

                    // Generate a vector of random cluster sizes
                    treatment_cluster_size = RandomClusterSize(sample_size[k], treatment_cluster_cum_vector);

                } else {

                    treatment_cluster_size = RandomClusterSize(sample_size[k], ExtractRow(treatment_cluster_cum_matrix, k - 1));

                }

                cluster_size_results.insert(cluster_size_results.end(), treatment_cluster_size.begin(), treatment_cluster_size.end());

            }

            treatment_size = treatment_cluster_size.size();

            // Number of patients in each cluster
            id_freq.clear();
            id_freq.insert(id_freq.end(), id_freq_control.begin(), id_freq_control.end());
            id_freq.insert(id_freq.end(), treatment_cluster_size.begin(), treatment_cluster_size.end());

            // Generate the data
            x.clear();
            y.clear();
            id.clear();

            x.insert(x.end(), x_control.begin(), x_control.end());
            y.insert(y.end(), y_control.begin(), y_control.end());
            id.insert(id.end(), id_control.begin(), id_control.end());

            // Current treatment arm
            for (i = 0; i < treatment_size; i++) {

                // Normal endpoint
                if (endpoint_index == 1) {

                    cluster_term = Normal(1, 0.0, between_cluster_sds[k])[0];
                    current_treatment_sample = Normal(treatment_cluster_size[i], means[k] + cluster_term, within_cluster_sds[k]);
                }

                // Binary endpoint
                if (endpoint_index == 2) {

                    cluster_term = Beta(1, treatment_alpha[k - 1], treatment_beta[k - 1])[0];
                    current_treatment_sample = Binary(treatment_cluster_size[i], cluster_term);

                }

                temp_vec = fillvec(treatment_cluster_size[i], 1.0);
                x.insert(x.end(), temp_vec.begin(), temp_vec.end());

                temp_vec = current_treatment_sample;
                y.insert(y.end(), temp_vec.begin(), temp_vec.end());

                temp_vec_int = FillVecInt(treatment_cluster_size[i], control_size + i + 1);
                id.insert(id.end(), temp_vec_int.begin(), temp_vec_int.end());

            }

            NumericMatrix x_mat(x.size(), 2);

            for (i = 0; i < x.size(); i++) {

                x_mat(i, 0) = 1.0;
                x_mat(i, 1) = x[i];

                // Direction of favorable outcome
                if (direction_index == 2) {
                    if (endpoint_index == 1) y[i] = -y[i];
                    if (endpoint_index == 2) y[i] = 1.0 - y[i];
                }

            }

            // Analyze the data for the current treatment-control comparison using generalized estimating equations

            // Normal endpoint
            if (endpoint_index == 1) results = ContGEE(y, x_mat, id, id_freq, direction_index, max_iterations, min_increment);

            // Binary endpoint
            if (endpoint_index == 2) results = BinGEE(y, x_mat, id, id_freq, direction_index, max_iterations, min_increment);

            // Save the convergence flags and raw p-values
            for (j = 0; j < 4; j++) pval(k - 1, j) = results.pval[j];

            // Save the convergence flags and effect estimates
            for (j = 0; j < 3; j++) eff(k - 1, j) = results.eff[j];

        }

        // Apply a multiplicity adjustment
        if (mult_test_index >= 1) {

            // Extract the raw pvalues for each method
            for (j = 1; j < 4; j++) {

                for (k = 0; k < narms - 1; k++) pvalue[k] = pval(k, j);

                adj_pvalue = TradMultAdj(mult_test_index, pvalue, weight, transition);

                for (k = 0; k < narms - 1; k++) pval(k, j) = adj_pvalue[k];

            }

        }

        // Save the convergence flags and adjusted p-values for each simulation run
        i = 0;
        for (k = 0; k < narms - 1; k++) {
            for (j = 0; j < 4; j++) {
                pval_results(sim, i) = pval(k, j);
                i++;
            }
        }                

        // Save the convergence flags and effect estimates for each simulation run
        i = 0;
        for (k = 0; k < narms - 1; k++) {
            for (j = 0; j < 3; j++) {
                coef_results(sim, i) = eff(k, j);
                i++;
            }
        }                

    }
    // End of simulations

    /*******************************************************************/

    // Return simulation results

    return List::create(Named("pval_results") = pval_results,
                        Named("coef_results") = coef_results,
                        Named("cluster_size_results") = cluster_size_results);


}
// End of ClustRandGEEC

// [[Rcpp::export]]
NumericVector ExportTradMultAdj(const int &test, const NumericVector &pvalue, const NumericVector &weight, const NumericVector &transition) {

    // Compute adjusted p-values
    vector<double> temp = TradMultAdj(test, FromNumericVector(pvalue), FromNumericVector(weight), FromNumericVector(transition));

    NumericVector adj_pvalue = ToNumericVector(temp);

    return adj_pvalue;

}

// [[Rcpp::export]]
IntegerVector ExportRandomClusterSize(const int &sample_size, const NumericVector &proportion) {

    vector<int> temp = RandomClusterSize(sample_size, FromNumericVector(proportion));

    IntegerVector cluster_size = ToIntegerVector(temp);

    return(cluster_size);

}
