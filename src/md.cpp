#include <Rcpp.h>
#include <RcppNumerical.h>

using namespace Numer;
using namespace std;
using namespace Rcpp;

#include "medstruct.h"
#include "medsupport.h"
#include "meddist.h"
#include "medstattest.h"
//#include "meddebug.h"

#include <exception>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <math.h> 
#include <iomanip>
#include <fstream>
#include <string>

int n_models = 4, direction_index_local = 1;
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

            if (pvalue[0] > alpha && pvalue[1] <= alpha / 2) { 
                outcome[0] = 0;
                outcome[1] = 1;
            }

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

        // Both populations are always selected if the interaction condition is not specified
        outcome3 = 1.0;     // # nocov

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

//REM, 43
double Logit(const double &x) {
    return log(x /(1.0 - x));
}

//REM, 47
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

                // Treatment selection rule
                if (futility_flag[i] < 1.0 && cp2[i] > max_cp) {
                    max_cp = cp2[i];
                    select_flag = i;    
                }

            }

            // Final analysis

            if (select_flag > -1) {

                control_sample = ExtractSamples(overall_control_sample, 0, sample_size_fa[0]);

                treatment_sample.clear();
                for (j = 0; j < sample_size_fa[select_flag + 1]; j++) treatment_sample.push_back(overall_treatment_sample(select_flag, j));

                if (endpoint_index == 1) test_result = TTest(control_sample, treatment_sample, 0.0, direction_index);
                if (endpoint_index == 2) test_result = PropTest(control_sample, treatment_sample, 0.0, direction_index);

                ad_outcome = (test_result.pvalue <= alpha / (narms - 1.0));        

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
                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

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
    
                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

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
    
                test_result = LogrankTest(outcome_censor_control, outcome_censor_treatment, 1.0, direction_index);

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
            }

            denom = 0.0;
            for (i = 0; i < n_doses - 1; i++) denom += pow(post_prob[i], balance);

            // Update the randomization ratio
            randomization_ratio[0] = ratio_placebo;
            for (i = 1; i < n_doses; i++) {
                randomization_ratio[i] = (1.0 - ratio_placebo) * pow(post_prob[i - 1], balance) / denom;
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
