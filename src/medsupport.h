// Supportive functions

#ifndef MEDSUPPORT_H
#define MEDSUPPORT_H

// Vector of constant values
vector<double> fillvec(const int &n, const double &value) {

    std::vector<double> result(n); 
    for (int i = 0; i < n; i++) {
        result[i] = value;
    }
    return result;  
}

// Vector of constant values
vector<int> FillVecInt(const int &n, const int &value) {

    std::vector<int> result(n); 
    for (int i = 0; i < n; i++) {
        result[i] = value;
    }
    return result;  
}

// Vector of treatment indicators
vector<int> FillTreatmentIndicators(const vector<int> &n) { 

    int i, n_patients = std::accumulate(n.begin(), n.end(), 0), n_arms = n.size(), start;

    vector<int> result(n_patients); 
    start = 0;
    for (i = 0; i < n_arms; i++) {
        fill(result.begin() + start, result.begin() + start + n[i], i); 
        start += n[i];
    }
    return result;  
}

// Extract a row
vector<double> ExtractRow(const NumericMatrix &mat, const int &index) {

    int i, m = mat.ncol();
    vector<double> row(m);
    for(i = 0; i < m; i++) row[i] = mat(index, i);

    return row;     

}

// Extract a column
vector<double> ExtractColumn(const NumericMatrix &mat, const int &index) {

    int i, m = mat.nrow();
    vector<double> column(m);
    for(i = 0; i < m; i++) column[i] = mat(i, index);

    return column;     

}

double Sq(const double &x) {

    return x * x;

}

// # nocov start
double Sign(const double &x) {

    double res;

    if (x < 0.0) {
        res = -1.0;
    } else {
        res = 1.0;        
    }

    return res;

}
// # nocov end

double sum(const vector<double> &vec) {
    int i, m = vec.size();
    double sum = 0.0;
    for(i = 0; i < m; ++i) sum += vec[i];
    return sum;
}

double sumsq(const vector<double> &vec) {
    int i, m = vec.size();
    double sum = 0.0;
    for(i = 0; i < m; ++i) sum += Sq(vec[i]);
    return sum;
}

double SumVec(const vector<double> &vec) {
    int i, m = vec.size();
    double sum = 0.0;
    for(i = 0; i < m; ++i) sum += vec[i];
    return sum;
}

int SumVecInt(const vector<int> &vec) {
    int i, m = vec.size();
    int sum = 0.0;
    for(i = 0; i < m; ++i) sum += vec[i];
    return sum;
}

// # nocov start
vector<double> vecsum(const vector<double> &x, const vector<double> &y) {
    int i, m = x.size();
    vector<double> sum(m);
    for(i = 0; i < m; ++i) sum[i] = x[i] + y[i];
    return sum;
}

double scalprod(const vector<double> &x, const vector<double> &y) {
    int i, m = x.size();
    double sum = 0.0;
    for(i = 0; i < m; ++i) sum += x[i] * y[i];
    return sum;
}


vector<double> AddVec(const vector<double> &x, const vector<double> &y) {
    int i, m = x.size();
    vector<double> sum(m);
    for(i = 0; i < m; ++i) sum[i] = x[i] + y[i];
    return sum;
}

// Compute averages by dividing by the number of simulations 
vector<double> ComputeAverage(vector<double> &vec, const int &nsims) {

    int i;
    int m = vec.size();
    vector<double> ave(m); 

    for (i = 0; i < m; i++) {
        ave[i] = vec[i] / nsims;
    }

    return ave;  

}
// # nocov end

#endif // MEDSUPPORT_H
