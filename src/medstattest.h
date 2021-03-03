// Statistical and related tests

#ifndef MEDSTATTEST_H
#define MEDSTATTEST_H


struct LogRankdata {
    double t;
    char cens;
    char id;
    bool operator < (const LogRankdata& str) const
    {
        return (t < str.t);
    }
    bool operator > (const LogRankdata& str) const
    {
        return (t >= str.t);
    }
};

TestResult PropTest(const vector<double> &x, const vector<double> &y, const double &margin, const int &direction) {

    double prop1 = 0.0, prop2 = 0.0;
    int i, n1 = x.size(), n2 = y.size();

    for (i = 0; i < n1; i++) prop1 += x[i];
    for (i = 0; i < n2; i++) prop2 += y[i];

    prop1 /= (n1 + 0.0);
    prop2 /= (n2 + 0.0);

    double ave = (prop1 * n1 + prop2 * n2)/ (n1 + n2 + 0.0);

    double se = sqrt(ave * (1 - ave) * ((1/(n1 + 0.0)) + (1/(n2 + 0.0))));

    double stat = (prop2 - prop1 + margin) / se;

    if (direction == 2) stat = - stat;

    if (std::isnan(stat)) {
        stat = -3.0;
    }

    TestResult res;

    // One-sided p-value
    res.test_stat = stat;
    res.pvalue = 1.0 - rcpp_pnorm(stat);

    return res;

}

TestResult TTest(const vector<double> &x, const vector<double> &y, const double &margin, const int &direction) {

    double Sm1=0;       // Sm1 = Sample 1 Mean.
    double Var1=0;       // Var1 = Sample 1  Variation.
    unsigned Sn1 = x.size();      // Sn1 = Sample 1 Size.
    double Sm2=0;       // Sm2 = Sample 2 Mean.
    double Var2=0;       // Var2 = Sample 2 Variation.
    unsigned Sn2 = y.size();     // Sn2 = Sample 2 Size.
    double j;
    for (int i = 0; i<Sn1 ; ++i)
    {
        j = x[i];
        Sm1 += j;
        Var1 += j*j;
    }
    Sm1 /=Sn1;
    Var1 = (Var1/Sn1-Sm1*Sm1)*Sn1/(Sn1-1);
    for (int i = 0; i<Sn2 ; ++i)
    {
        j = y[i];
        Sm2 += j;
        Var2 += j*j;
    }
    Sm2 /=Sn2;
    Sm2 -= margin;
    Var2 = (Var2/Sn2-Sm2*Sm2)*Sn2/(Sn2-1);
    // Degrees of freedom:
    double v = Sn1 + Sn2 - 2;
    // Pooled variance:
    double sp = sqrt(((Sn1-1) * Var1  + (Sn2-1) * Var2 ) / v);
    // t-statistic:
    double stat = (Sm2 - Sm1) / (sp * sqrt(1.0 / Sn1 + 1.0 / Sn2));
    if (direction == 2) stat = - stat;

    TestResult res;

    if (std::isnan(stat)) {
        stat = -3.0;
    }

    // One-sided p-value based on a t distribution
    res.test_stat = stat;
    res.pvalue = 1.0 - rcpp_pt(stat, v);

    return res;

}


void TupleSort(const std::vector<double> &in1, const std::vector<char> &in2, const std::vector<char> &in3, vector<LogRankdata> &vec) {

    vec.resize(in1.size());
    for (int i = 0; i < in1.size();++i){
        vec[i].t = in1[i];
        vec[i].cens = in2[i];
        vec[i].id = in3[i];
    }
    sort(vec.begin(),vec.end());
}


TestResult CoreLogrankTest(const std::vector<double> &xv, const std::vector<double> &yv, const std::vector<char> &cxv, const std::vector<char> &cyv, const double &margin, const int &direction) {

    int i;

    std::vector<char> ixv(xv.size(),1);
    std::vector<char> iyv(yv.size(),2);

    std::vector<double> xy;
    xy.reserve(xv.size()+yv.size());
    xy.insert(xy.end(),xv.begin(),xv.end());
    xy.insert(xy.end(),yv.begin(),yv.end());

    std::vector<char> cxy;
    cxy.reserve(cxv.size()+cyv.size());
    cxy.insert(cxy.end(),cxv.begin(),cxv.end());
    cxy.insert(cxy.end(),cyv.begin(),cyv.end());

    std::vector<char> ixy;
    ixy.reserve(ixv.size()+iyv.size());
    ixy.insert(ixy.end(),ixv.begin(),ixv.end());
    ixy.insert(ixy.end(),iyv.begin(),iyv.end());

    std::vector<LogRankdata> vec;

    TupleSort(xy,cxy,ixy,vec);

    std::vector<double> t;
    t.reserve(vec.size());
    std::vector<int> m1;
    m1.reserve(vec.size());
    std::vector<int> m2;
    m2.reserve(vec.size());
    std::vector<int> cens1,cens2;
    cens1.reserve(vec.size());cens2.reserve(vec.size());
    int curit=0;
    for (i =0; i <vec.size();++i)
    {
        t.push_back(vec[i].t);
        m1.push_back(0);
        m2.push_back(0);
        cens1.push_back(0);
        cens2.push_back(0);
        if(vec[i].id == 1)
        {
            m1[curit] += 1-vec[i].cens;
            cens1[curit]+=vec[i].cens;
        }
        else
        {
            m2[curit] += 1-vec[i].cens;
            cens2[curit]+=vec[i].cens;
        }
        while (vec[i].t==vec[i+1].t)
        {
            ++i;
            if(vec[i].id == 1)
            {
                m1[curit] += 1-vec[i].cens;
                cens1[curit]+=vec[i].cens;
            }
            else
            {
                m2[curit] += 1-vec[i].cens;
                cens2[curit] += (int)vec[i].cens;
            }
        }
        ++curit;
    }

    int n1 = xv.size(), n2 = yv.size();
    double s1 = 0, s2 = 0, oe = 0, v = 0;


    for(int i = 0; i<curit;++i)
    {
        oe = m1[i]-((double)n1)/(n1+(n2+0.0)*margin);
        v = (double)n1*(n2+0.0)*margin/((n1+(n2+0.0)*margin)*(n1+(n2+0.0)*margin));
        v = (v!=v)?0:v;
        s1 += (m1[i]+m2[i])*oe;
        s2 += (m1[i]+m2[i])*v;
        n1 -= m1[i]+cens1[i];
        n2 -= m2[i]+cens2[i];
    }

    double stat = s1/sqrt(s2);

    if (std::isnan(stat)) {
        stat = -5.0;
    }
    if (direction == 2) stat = - stat;

    TestResult res;

    // One-sided logrank p-value
    res.test_stat = stat;
    res.pvalue = 1.0 - rcpp_pnorm(stat);

    return res;

}

double EventCount(const OutcomeCensor &outcome_censor_x, const OutcomeCensor &outcome_censor_y) {
    
    // Extract outcome and censor objects
    vector<double> x = outcome_censor_x.outcome;
    vector<double> cx = outcome_censor_x.censor;

    vector<double> y = outcome_censor_y.outcome;
    vector<double> cy = outcome_censor_y.censor;

    double count = cx.size() + cy.size() + 0.0 - std::accumulate(cx.begin(), cx.end(), 0.0) - std::accumulate(cy.begin(), cy.end(), 0.0); 

    return count;
}


TestResult LogrankTest(const OutcomeCensor &outcome_censor_x, const OutcomeCensor &outcome_censor_y, const double &margin, const int &direction) {

    // Extract outcome and censor objects
    vector<double> x = outcome_censor_x.outcome;
    vector<double> cx = outcome_censor_x.censor;

    vector<double> y = outcome_censor_y.outcome;
    vector<double> cy = outcome_censor_y.censor;

    vector<char> cx1(cx.size());
    std::transform(cx.begin(),cx.end(),cx1.begin(),[](double d){return d>0;});

    vector<char> cy1(cy.size());
    std::transform(cy.begin(),cy.end(),cy1.begin(),[](double d){return d>0;});

    TestResult res;

    res = CoreLogrankTest(x, y, cx1, cy1, margin, direction);

    return res;
}

double FindMilestone(const vector<int> &stratum_list, const vector<int> &stratum, const vector<double> &local_start, const int &target) {

    double milestone = 0.0;
    vector<double> vec;
    int i, m = local_start.size();

    // Remove the times corresponding to censored observations
    for(i = 0; i < m; i++) {
        if (local_start[i] >= 0.0 && find(stratum_list.begin(), stratum_list.end(), stratum[i]) != stratum_list.end()) vec.push_back(local_start[i]);
    } 

    if (vec.size() > 0) {

        sort(vec.begin(), vec.end());

        if (target > vec.size()) {
            milestone = vec.back();
        } else {        
            milestone = vec[target - 1];
        }

    } else {
        milestone = 10000.0;        
    }


    return milestone;    

}

// Extract subvectors of outcomes and censoring indicators from vectors of outcomes and censoring indicators
OutcomeCensor ExtractOutcomeCensor(const vector<int> &stratum_list, const vector<int> &stratum, const vector<double> &start, const vector<double> &local, const vector<double> &local_censor, const double &milestone) {

    OutcomeCensor outcome_censor;
    vector<double> outcome, censor;
    int i, m = start.size();
    double temp, temp_censor;

    for (i = 0; i < m; i++) { 

        if (start[i] <= milestone) {

            temp_censor = local_censor[i];
            if (start[i] + local[i] <= milestone) {
                temp = local[i];
            } else {
                temp = milestone - start[i];
                temp_censor = 1.0;
            }
            if (find(stratum_list.begin(), stratum_list.end(), stratum[i]) != stratum_list.end()) {
                outcome.push_back(temp);
                censor.push_back(temp_censor);
            }
        }

    }

    outcome_censor.outcome = outcome;
    outcome_censor.censor = censor;

    return outcome_censor;

}

// Estimate the hazard ratio assuming an exponential distribution
double HazardRatio(const OutcomeCensor &outcome_censor_control, const OutcomeCensor &outcome_censor_treatment) {

    double rate_control, rate_treatment, hazard_ratio;    
    int m_control, m_treatment;

    // Hazard rates
    m_control = outcome_censor_control.outcome.size();
    rate_control = (m_control + 0.0 - sum(outcome_censor_control.censor)) / sum(outcome_censor_control.outcome);

    m_treatment = outcome_censor_treatment.outcome.size();
    rate_treatment = (m_treatment + 0.0 - sum(outcome_censor_treatment.censor)) / sum(outcome_censor_treatment.outcome);

    // Hazard ratio
    hazard_ratio = rate_treatment / rate_control;

    return hazard_ratio;  

}

// Combination function based on stagewise test statistics
double CombFunctionTestStat(const double &test_stat1, const double &test_stat2, const double &w1, const double &w2) {

    double pvalue = 1.0 - rcpp_pnorm(sqrt(w1) * test_stat1 + sqrt(w2) * test_stat2);

    return pvalue;

}

vector<double> MarginalCombTest(const double &p1, const double &p2, const int &test) {

    vector<double> adjusted_pvalue(2);
    double intersection_pvalue = 1.0;

    // Holm test
    if (test == 1) intersection_pvalue = 2.0 * min(p1, p2); 

    adjusted_pvalue[0] = max(p1, intersection_pvalue);
    adjusted_pvalue[1] = max(p2, intersection_pvalue);

    return adjusted_pvalue;

}

double IntersectionPvalue(const double &p1, const double &p2, const int &test) {

    double intersection_pvalue = 1.0;

    // Holm test
    if (test == 1) intersection_pvalue = 2.0 * min(p1, p2); 

    // Hochberg test
    if (test == 2) intersection_pvalue = 2.0 * min(p1, p2); 

    return intersection_pvalue;

}


#endif // MEDSTATTEST_H
