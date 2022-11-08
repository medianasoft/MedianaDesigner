// Distribution functions and generation of random variables

#ifndef MEDDIST_H
#define MEDDIST_H

double rcpp_pnorm(const double &x) {
    NumericVector vec_input(1), vec_output(1);
    vec_input[0] = x;
    vec_output = Rcpp::pnorm(vec_input);
    return vec_output[0];
}

double rcpp_qnorm(const double &x) {
    NumericVector vec_input(1), vec_output(1);
    vec_input[0] = x;
    vec_output = Rcpp::qnorm(vec_input);
    return vec_output[0];
}

double rcpp_pt(const double &x, const double &df) {
    NumericVector vec_input(1), vec_output(1);
    vec_input[0] = x;
    vec_output = Rcpp::pt(vec_input, df);
    return vec_output[0];
}

// Truncated exponential CDF 
double TruncatedExponentialCDF(const double &x, const double &par) {

    double y;

    if (abs(par) < 0.0001) y = x; else y = (1.0 - exp(- par * x)) / (1.0 - exp(- par));

    return y;

}

// Find the parameter for a given median enrollment time 
double TruncatedExponentialMedian(const double &median) {

        double lower_bound, midpoint, upper_bound;

        lower_bound = -9.0;
        upper_bound = 10.0;
        midpoint = (lower_bound + upper_bound) / 2.0;

        do {

            if (TruncatedExponentialCDF(median, midpoint) < 0.5) 
                { 
                    lower_bound = midpoint;
                } else {   
                    upper_bound = midpoint;        
                }
            midpoint = (lower_bound + upper_bound) / 2.0;

        }
        while (upper_bound - lower_bound >= 0.001);

        return midpoint;

}


// Vector of exponentially distributed values 
vector<double> Exponential(const int &n, const double &lambda) {

    NumericVector temp_vector = Rcpp::rexp(n, lambda);
    vector<double> result = as<vector<double>>(temp_vector); 
    return result;  

}

// Vector of binary values 
vector<double> Binary(const int &n, const double &prop) {

    NumericVector result = Rcpp::rbinom(n, 1, prop);
    vector<double> result_vector = as<vector<double>>(result);
    return result_vector;  

}

// Vector of binary values 
vector<double> Normal(const int &n, const double &mean, const double &sd) {

    NumericVector result = Rcpp::rnorm(n, mean, sd);
    vector<double> result_vector = as<vector<double>>(result);
    return result_vector;  

}

// # nocov start
// Vector of uniformly distributed values 
vector<double> Uniform(const int &n, const double &min, const double &max) {

    NumericVector temp_vector = Rcpp::runif(n, min, max);
    vector<double> result = as<vector<double>>(temp_vector); 
    return result;  

}
// # nocov end

// Vector of gamma distributed values 
vector<double> Gamma(const int &n, const double &shape, const double &rate) {

    NumericVector temp_vector = Rcpp::rgamma(n, shape, 1.0 / rate);
    vector<double> result = as<vector<double>>(temp_vector); 
    return result;  

}

// Vector of beta distributed values 
vector<double> Beta(const int &n, const double &shape1, const double &shape2) {

    NumericVector temp_vector = Rcpp::rbeta(n, shape1, shape2);
    vector<double> result = as<vector<double>>(temp_vector); 
    return result;  

}

// Vector of truncated exponential values 
vector<double> TruncatedExponential(const int &n, const double &par, const double &min, const double &max) {

    NumericVector temp_vector(n);
    double temp;
    int i;

    if (par == 0.0) {
        temp_vector = Rcpp::runif(n, min, max);     // # nocov
    } else {
        for (i = 0; i < n; i++) {
            temp = - log(1.0 - Rcpp::runif(1, 0.0, 1.0)[0] * (1.0 - exp(-par))) / par;
            temp_vector[i] = min + temp * (max - min);
        }
    }

    vector<double> result = as<vector<double>>(temp_vector); 

    return result;  

}

// Vector of enrollment times 
vector<double> Enrollment(const int &n, const double &enrollment_period, const int &enrollment_distribution, const double &enrollment_parameter) {

    vector<double> result(n); 
    double trunc_exp_par = 0.0;

    // Uniform enrollment
    if (enrollment_distribution == 1) { 

        result = Uniform(n, 0.0, enrollment_period);    // # nocov

    }

    // Truncated exponential enrollment distribution
    if (enrollment_distribution == 2) {

        trunc_exp_par = TruncatedExponentialMedian(enrollment_parameter / enrollment_period);
        result = TruncatedExponential(n, trunc_exp_par, 0.0, enrollment_period);

    }


    return result;  

}

// Vector of times to dropout
vector<double> Dropout(const int &n, const int &dropout_distribution, const vector<double> &dropout_parameter) {

    vector<double> result(n); 

    double dropout_rate;

    // No dropout
    if (dropout_distribution == 1) { 

        result = fillvec(n, 100000.0);      // # nocov

    }

    // Exponential dropout
    if (dropout_distribution == 2) { 

        dropout_rate = -log(1.0 - dropout_parameter[0])/dropout_parameter[1];
        result = Exponential(n, dropout_rate);

    }

    return result;  

}

// # nocov start
// Multivariate normal distribution (single vector)
vector<double> MVNormal(const int &m, const vector<double> &mean, const vector<double> &sd, const NumericMatrix &corr) {

    int i, j, k;
    double csum;

    vector<double> normal_data(m), mv_normal_data(m);
    NumericMatrix chol(m, m);

    for(i = 0; i< m; i++) {
        for(j = 0; j < m; j++) {
            chol(i, j) = 0.0;
        }
    }    

    // Cholesky lower diagonal matrix
    // http://rosettacode.org/wiki/Cholesky_decomposition
    for(i = 0; i< m; i++) {
        for(j = 0; j < (i + 1); j++) {
            csum = 0.0;
            for(k = 0; k < j; k++){
                csum += chol(i, k) * chol(j, k);
            }
            if (i == j) {
                chol(i, j) = sqrt(corr(i, i) - csum);
            } else {
                chol(i, j) = (corr(i, j) - csum) / chol(j, j);
            }
        }
    }

    // Generate uncorrelated normal vectors
    for (i = 0; i< m; i++) {
        normal_data[i] = Normal(1, 0.0, 1.0)[0];
    }

    // Compute correlated normal vectors
    for(i = 0; i< m; i++) {
            csum = 0;
            for(k = 0; k < m; k++) {
                csum += chol(i, k) * normal_data[k];
            }
            mv_normal_data[i] = csum * sd[i] + mean[i];
    }

    return mv_normal_data;

}

vector<double> MVNormalRho(const int &m, const vector<double> &mean, const vector<double> &sd, const double &rho) {

    NumericMatrix corr(m, m);

    for(unsigned i = 0; i < m; i++){
        for(unsigned j = 0; j < m; j++){
            if (i == j) corr(i, j) = 1.0; else corr(i, j) = rho;
        }
    }

    return MVNormal(m, mean, sd, corr);
}
// # nocov end

#endif // MEDDIST_H
