// Supportive functions
// # nocov start

#ifndef MEDDEBUG_H
#define MEDDEBUG_H

void printVectorDouble(const string &label, const vector<double> &vec) {

    int m = vec.size();

    cerr<<label<<": ";
    for(int i = 0; i < m; ++i) {
        cerr<<vec[i]<<", ";
    }
    cerr<<endl;

}

void printVectorInt(const string &label, const vector<int> &vec) {

    int m = vec.size();

    cerr<<label<<": ";

    if (m > 0) {
        for(int i = 0; i < m; ++i) {
            cerr<<vec[i]<<", ";
        }
    } else {
        cerr<<"NA";
    }
    cerr<<endl;

}

void printMatrixDouble(const string &label, const NumericMatrix &mat) {

    int i, j, m1 = mat.nrow(), m2 = mat.ncol();

    cerr<<label<<": "<<endl;
    for(i = 0; i < m1; i++) {
       for(j = 0; j < m2; j++) {
            cerr<<mat(i, j)<<", ";
        }
        cerr<<endl;
    }
    cerr<<endl;

}

#endif // MEDDEBUG_H
// # nocov end