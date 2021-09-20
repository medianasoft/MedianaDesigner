#ifndef MEDDESMULT_H
#define MEDDESMULT_H

typedef std::pair<double,double> weighted_pvalue;

// # nocov start
bool SortByFirstWeighted (const weighted_pvalue& l, const weighted_pvalue& r) { 
    return l.first < r.first; 
}

bool SortBySecondWeighted (const weighted_pvalue& l, const weighted_pvalue& r) { 
    return l.second < r.second; 
}
// # nocov end

int ArgMin(const std::vector<double> &p, const std::vector<double> &w, const std::vector<int> &processed) {

  int i, index = -1, m = p.size();
  double pmin;

  for(i = 0; i < m; i++) {
    if (w[i] > 0.0 && processed[i] == 0) {
      if (index == -1) {
        pmin = p[i]/w[i];
        index = i;
      }
      if (index >= 0 && p[i]/w[i] < pmin && processed[i] == 0) {
        pmin = p[i]/w[i];
        index = i;
      }
    }
  }
  return index;
}

vector<double> ChainAdj(const std::vector<double> &p, const std::vector<double> &par1, const std::vector<double> &par2) {

    int m = p.size();
    int i, j, k, ind;
    vector<double> w(m), w_temp(m), adjpvalue(m);
    vector<int> processed(m);
    NumericMatrix g(m, m), g_temp(m, m);
    double pmax = 0.0; 

    for(i = 0; i < m; i++) w[i] = par1[i];

    k = 0;    
    for(i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {
            g(i, j) = par2[k];
            k++;
        }
    }

    for(i = 0; i < m; i++) {
      adjpvalue[i] = 0.0;  
      processed[i] = 0;         
    }

    // Loop over all null hypotheses
    for(i = 0; i < m; i++) {
      // Find the index of the smallest weighted p-value among the non-processed null hypotheses
      ind = ArgMin(p, w, processed);

      if (ind>-1) {
        adjpvalue[ind] = max(p[ind]/w[ind], pmax);
        adjpvalue[ind] = min(1.0, adjpvalue[ind]);
        pmax = adjpvalue[ind];
        // This null hypothesis has been processed
        processed[ind] = 1;

        // Update the hypothesis weights after a null hypothesis has been processed
        for(j = 0; j < m; j++) w_temp[j] = w[j];
        for(j = 0; j < m; j++) {
          if (processed[j] == 0)
            w[j] = w_temp[j] + w_temp[ind] * g(ind, j); else w[j] = 0.0;
        }

        // Update the transition parameters (connection weights) after the rejection
        for(j = 0; j < m; j++) {
          for(k = 0; k < m; k++) {
                    g_temp(j, k) = g(j, k);
            }
        }

        for(j = 0; j < m; j++) {
            for(k = 0; k < m; k++) {
                if (processed[j] == 0 && processed[k] == 0 && j != k && g_temp(j, ind) * g_temp(ind, j) != 1.0)
                    g(j, k) = (g_temp(j, k) + g_temp(j, ind) * g_temp(ind, k))/(1 - g_temp(j, ind) * g_temp(ind, j)); else g(j, k) = 0.0;
                }
            }
        } else {
            // # nocov start
            for(j = 0; j < m; j++) {
                if (processed[j] == 0) adjpvalue[j] = 1.0;
            }
            // # nocov end
        }
    }

    return adjpvalue;

}

// # nocov start
std::vector<double> PartialSum(const std::vector<double> &x) {

    int i, m = x.size();
    vector<double> partial(m);

    if (m == 1) partial[0] = x[0];

    if (m > 1) {
        
        partial[0] = x[0];
        for(i = 1; i < m; i++) {

            partial[i] = x[i] + partial[i - 1];

        }
    }

    return partial;

}
// # nocov end

double BonferroniGlobal(const std::vector<double> &pvalue, const std::vector<double> &weight) {

    int i, m = pvalue.size(), nonzero_m;
    double globalp = 1.0;
    vector<double> nonzero_p, nonzero_w;

    // Remove zero weights and count the number of non-zero weights 
    nonzero_m = 0;
    nonzero_p.clear();
    nonzero_w.clear();
    for(i = 0; i < m; i++) {

        if (abs(weight[i]) > 0.000001) {
            nonzero_p.push_back(pvalue[i]);
            nonzero_w.push_back(weight[i]);
            nonzero_m++;
        }

    }

    if (nonzero_m > 0) {

        globalp = nonzero_p[0] / nonzero_w[0];

        for(i = 1; i < nonzero_m; i++) {
            globalp = min(globalp, nonzero_p[i] / nonzero_w[i]);
        }    

    }

    return globalp;


}

// # nocov start
double SimesGlobal(const std::vector<double> &pvalue, const std::vector<double> &weight) {
    
    int i, m = pvalue.size(), nonzero_m;
    double globalp = 1.0;
    vector<double> nonzero_p, nonzero_w, cumsum;
    vector<weighted_pvalue> ordered_p;
    
    // Remove zero weights and count the number of non-zero weights 
    nonzero_m = 0;
    nonzero_p.clear();
    nonzero_w.clear();
    for(i = 0; i < m; i++) {

        if (abs(weight[i]) > 0.000001) {
            nonzero_p.push_back(pvalue[i]);
            nonzero_w.push_back(weight[i]);
            nonzero_m++;
        }

    }

    if (nonzero_m > 1) {

        vector<double> modu(nonzero_m), temp_vec(nonzero_m), cumsum(nonzero_m);

        // Sort by the p-value
        ordered_p.clear();
        for(i = 0; i < nonzero_m; i++) {
            ordered_p.push_back(weighted_pvalue(nonzero_p[i], nonzero_w[i]));
        }    

        // Sort by p-value
        sort(ordered_p.begin(), ordered_p.end(), SortByFirstWeighted);

        for(i = 0; i < nonzero_m; i++) {
            modu[i] = ordered_p[i].second;
        }    

        cumsum = PartialSum(modu);

        for(i = 0; i < nonzero_m; i++) {
            temp_vec[i] = ordered_p[i].first / cumsum[i];
        }    

        globalp = *std::min_element(temp_vec.begin(), temp_vec.end());    

    }

    if (nonzero_m == 1) {
        globalp = min(1.0, nonzero_p[0] / nonzero_w[0]);
    }


    return globalp;    

}

double IncompleteSimesGlobal(const std::vector<double> &pvalue, const std::vector<double> &weight) {
    
    int i, m = pvalue.size(), nonzero_m;
    double globalp = 1.0;
    vector<double> nonzero_p, nonzero_w, cumsum;
    vector<weighted_pvalue> ordered_p;
    
    // Remove zero weights and count the number of non-zero weights 
    nonzero_m = 0;
    nonzero_p.clear();
    nonzero_w.clear();
    for(i = 0; i < m; i++) {

        if (abs(weight[i]) > 0.000001) {
            nonzero_p.push_back(pvalue[i]);
            nonzero_w.push_back(weight[i]);
            nonzero_m++;
        }

    }

    if (nonzero_m > 1) {

        vector<double> modu(nonzero_m), temp_vec(nonzero_m), cumsum(nonzero_m);

        // Sort by the p-value
        ordered_p.clear();
        for(i = 0; i < nonzero_m; i++) {
            ordered_p.push_back(weighted_pvalue(nonzero_p[i], nonzero_w[i]));
        }    

        // Sort by p-value
        sort(ordered_p.begin(), ordered_p.end(), SortByFirstWeighted);

        modu[0] = 0.0;

        for(i = 1; i < nonzero_m; i++) {
            modu[i] = ordered_p[i - 1].second;
        }    

        cumsum = PartialSum(modu);

        for(i = 0; i < nonzero_m; i++) {
            temp_vec[i] = (1.0 - cumsum[i]) * ordered_p[i].first / ordered_p[i].second;
        }    

        globalp = *std::min_element(temp_vec.begin(), temp_vec.end());    

    }

    if (nonzero_m == 1) {
        globalp = min(1.0, nonzero_p[0] / nonzero_w[0]);
    }


    return globalp;    

}
// # nocov end

// Weighted Holm, Hochberg and Hommel tests implemented using the closed testing principle
vector<double> ClosedTestingAdj(const int &test, const std::vector<double> &pvalue, const std::vector<double> &weight) {

    int i, j, k, m = pvalue.size();
    vector<double> adjpvalue(m), loc_weight(m);
    int nint = (int) pow (2, m) - 1;
    double weight_sum, globalp;

    NumericMatrix intersection(nint, m), intp(nint, m);

    // Populate the intersection matrix
    for (j = 0; j < nint; ++j) {
        
        for (i = 0; i < m; ++i) {
            k = floor(j/pow(2, m - i - 1));
            intersection(j, i) = 0.0;
            if (k/2.0 == floor(k/2.0)) {
                intersection(j, i) = 1.0;
            }
        }

        weight_sum = 0.0; 
        for (i = 0; i < m; i++) weight_sum += (intersection(j, i) * weight[i]);

        if (weight_sum > 0) {    

            for (i = 0; i < m; i++) loc_weight[i] = intersection(j, i) * weight[i] / weight_sum;

            // Holm
            if (test == 2) globalp = BonferroniGlobal(pvalue, loc_weight);     
            // Hochberg
            if (test == 3) globalp = IncompleteSimesGlobal(pvalue, loc_weight);     
            // Hommel
            if (test == 4) globalp = SimesGlobal(pvalue, loc_weight);     

        } else {

            globalp = 1.0;  // # nocov

        }

        for (i = 0; i < m; i++) intp(j, i) = globalp * intersection(j, i);

    }

    // Compute adjusted p-values
    for (i = 0; i < m; i++) {
        adjpvalue[i] = 0.0;
        for (j = 0; j < nint; j++) {
            adjpvalue[i] = max(adjpvalue[i], intp(j, i));
        }
    }   

    return adjpvalue;

}

// # nocov start
vector<double> BonferroniAdj(const std::vector<double> &pvalue, const std::vector<double> &weight) {

    int i, m = pvalue.size();
    vector<double> adjpvalue(m);

    for(i = 0; i < m; i++) {
        if (weight[i] > 0.0) adjpvalue[i] = min(1.0, pvalue[i] / weight[i]); else adjpvalue[i] = 1.0;
    }   

    return adjpvalue;

}
// End of BonferroniAdj
// # nocov end

vector<double> FixedSeqAdj(const std::vector<double> &pvalue, const std::vector<double> &weight) {

    int i, index, m = pvalue.size();
    vector<double> adjpvalue(m);
    vector<int> order(m);

    double current_adj_pvalue;

    // Define order
    for (i = 0; i < m; i++) order[i] = (int) (weight[i] - 1.0);

    for (i = 0; i < m; i++) {
        index = order[i]; 
        if (i == 0) {
            adjpvalue[index] = pvalue[index];
            current_adj_pvalue = adjpvalue[index];
        }
        if (i > 0) {
            adjpvalue[index] = max(pvalue[index], current_adj_pvalue);
            current_adj_pvalue = adjpvalue[index];
        }
    }    

    return adjpvalue;


}
// End of FixedSeqAdj

// Global test based on the truncated Hochberg procedure with equal weights
double HochbergGlobal(const std::vector<double> &pvalue, const int &n, const double &gamma) {

    int m = pvalue.size(), i;
    double globalp = 1.0;
    vector<double> p(pvalue);

    if (m > 0 && n > 0)  {

        // Truncated Hochberg procedure with gamma>0    
        if (gamma > 0.0) {

            vector<double> denom(m), sortp(m);

            sort(p.begin(), p.end());

            for(i = 0; i < m; i++) {
                sortp[i] = p[i] / (gamma / (m - (i + 1.0) + 1.0) + (1.0 - gamma) / (n + 0.0));    
            }

            globalp = *std::min_element(sortp.begin(), sortp.end());   
            
        }
        
        // Bonferroni test with gamma=0
        if (gamma == 0.0) {

            globalp = (n + 0.0) * (*std::min_element(p.begin(), p.end()));      // # nocov

        }

    }

    return globalp;

}
// End of HochbergGlobal

// Global test based on the truncated Hommel procedure with equal weights
double HommelGlobal(const std::vector<double> &pvalue, const int &n, const double &gamma) {

    int m = pvalue.size(), i;
    double globalp = 1.0;
    vector<double> p(pvalue);

    if (m > 0 && n > 0)  {

        // Truncated Hommel procedure with gamma>0    
        if (gamma > 0.0) {

            vector<double> denom(m), sortp(m);

            sort(p.begin(), p.end());

            for(i = 0; i < m; i++) {
                sortp[i] = p[i] / ((i + 1.0) * gamma / (m + 0.0) + (1.0 - gamma) / (n + 0.0));    
            }

            globalp = *std::min_element(sortp.begin(), sortp.end());   
            
        }
        
        // Bonferroni test with gamma=0
        if (gamma == 0.0) {

            globalp = (n + 0.0) * (*std::min_element(p.begin(), p.end()));      // # nocov

        }

    }

    return globalp;

}
// End of HommelGlobal

// Compute intersection p-value 
// PVALUE: Vector of component p-values.
// C: Vector of family weights.
double IntPvalue(const vector<double> &pvalue, const vector<double> &c) {

    double intp = 1.0;

    int m = pvalue.size(), i;

    for(i = 0; i < m; i++) {
      if (c[i] > 0.0 && pvalue[i]/c[i] < intp) intp = pvalue[i] / c[i];
    }

    return intp;

}


// Evaluate error fraction function for a family based on Hochberg or Hommel procedures 
// K: Number of null hypotheses included in the intersection within the family.
// N: Total number of null hypotheses in the family.  
// GAMMA: Truncation parameter (0<=GAMMA<1).
double ErrorFrac(const double &k, const double &n, const double &gamma) {
  
  double f = 0.0;

  if (k>0) {
    if (k == n) f = 1.0; else f = gamma + (1.0 - gamma) * (k + 0.0)/ (n + 0.0);
  }
  
  return f;
  
}

// Computation of adjusted p-values for gatekeeping procedures based on the standard and modified mixture methods
vector<double> MixtureProcAdjP(const int &nfam, const int &nperfam, const vector<double> &pvalue, const int &proc, const vector<double> &gamma, const int &method) {

    int nhyp = pvalue.size(), i, j, k, l, m, tot;
    int nint = (int) pow (2, nhyp) - 1;
    double sum;

    vector<double> adjpvalue(nhyp), temp_vec, pselected, intersection;

    NumericMatrix int_orig(nint, nhyp), int_rest(nint, nhyp), fam_rest(nint, nhyp), intp(nint, nhyp);

    // Construct the intersection index sets (int_orig) before the logical restrictions are applied. 
    // Each row is a vector of binary indicators (1 if the hypothesis is included in the original index set and 0 otherwise)       
    // Construct the intersection index sets (int_rest) and family index sets (fam_rest) after the logical restrictions are applied.   
    // Each row is a vector of binary indicators (1 if the hypothesis is included in the restricted index set and 0 otherwise)  
    for (j = 0; j < nint; ++j) {
        for (i = 0; i < nhyp; ++i) {
            k = floor(j/pow(2, nhyp - i - 1));
            int_orig(j, i) = 0.0;
            int_rest(j, i) = 0.0;
            if (k/2.0 == floor(k/2.0)) {
                int_orig(j, i) = 1.0;
                int_rest(j, i) = 1.0;
            }    
            fam_rest(j, i) = 1;
            intp(j, i) = 1;
        }
    }

    for (i = 0; i < nint; i++) {
    for (j = 0; j < nfam - 1; j++)  {
      for (k = 0; k < nperfam; k++)
      {
        // Index of the current null hypothesis in Family j
        m = j * nperfam + k + 1;
        // If this null hypothesis is included in the intersection hypothesis
        // all dependent null hypotheses must be removed from the intersection hypothesis
        if (int_orig(i, m - 1) == 1) {
          for (l = 0; l < nfam - (j + 1); l++) {
            int_rest(i, m + (l + 1) * nperfam - 1) = 0;        
            fam_rest(i, m + (l + 1) * nperfam - 1) = 0;                    
          }
        }
      }
    }
    }
  
    // Number of null hypotheses from each family included in each intersection 
    // before the logical restrictions are applied
    NumericMatrix korig(nint,nfam);
      
    // Number of null hypotheses from each family included in the current intersection
    // after the logical restrictions are applied
    NumericMatrix krest(nint,nfam); 
      
    // Number of null hypotheses from each family after the logical restrictions are applied
    NumericMatrix nrest(nint,nfam);   

    // Compute korig, krest and nrest
    for (i = 0; i < nint; i++) {    
        for (j = 0; j < nfam; j++) {

          // Index vector in the current family
          temp_vec = ExtractRow(int_orig, i);
          sum = 0.0;
          for (k = j * nperfam; k < (j + 1) * nperfam; k++) sum += temp_vec[k];
          korig(i,j) = sum;

          temp_vec = ExtractRow(int_rest, i);
          sum = 0.0;
          for (k = j * nperfam; k < (j + 1) * nperfam; k++) sum += temp_vec[k];
          krest(i,j) = sum;

          temp_vec = ExtractRow(fam_rest, i);
          sum = 0.0;
          for (k = j * nperfam; k < (j + 1) * nperfam; k++) sum += temp_vec[k];
          nrest(i,j) = sum;

        }
    } 

    // Index of the last non-empty restricted subset for each intersection
    vector<double> last_nonempty(nint);

    for (i = 0; i < nint; i++) {
        last_nonempty[i] = -1;    
        for (j = 0; j < nfam; j++) {
            if (krest(i,j) > 0) last_nonempty[i] = j;
        }
    }

    // Vector of intersection p-values
    vector<double> pint(nint);

    // Matrix of component p-values within each intersection
    NumericMatrix pcomp(nint,nfam);

    // Matrix of family weights within each intersection 
    NumericMatrix c(nint,nfam);

    // P-value for each hypothesis within each intersection    
    NumericMatrix p(nint,nhyp);

      // Compute the intersection p-value for each intersection hypothesis
      for (i = 0; i < nint; i++) {
        
        // Compute component p-values
        for(j = 0; j < nfam; j++) {

          // Consider non-empty restricted index sets 
          if (krest(i,j) > 0) {

            // Restricted index set in the current family 
            temp_vec = ExtractRow(int_rest, i);
            pselected.clear();
            for (k = j * nperfam; k < (j + 1) * nperfam; k++) {
                if (temp_vec[k] == 1) pselected.push_back(pvalue[k]);
            }

            // Total number of hypotheses used in the computation of the component p-value

            // Standard or modified methods
            if (method == 1 || method == 2) {
                if (method == 1) tot = nperfam;
                if (method == 2) tot = (int) nrest(i,j);
                if (proc == 1) pcomp(i,j) = HochbergGlobal(pselected, tot, gamma[j]);
                if (proc == 2) pcomp(i,j) = HommelGlobal(pselected, tot, gamma[j]);
            }

            // Enhanced method
            if (method == 3) {
                // # nocov start
                tot = (int) nrest(i,j);
                if (last_nonempty[i] != j) {
                    if (proc == 1) pcomp(i,j) = HochbergGlobal(pselected, tot, gamma[j]);
                    if (proc == 2) pcomp(i,j) = HommelGlobal(pselected, tot, gamma[j]);
                } else {
                    if (proc == 1) pcomp(i,j) = HochbergGlobal(pselected, tot, 1.0);
                    if (proc == 2) pcomp(i,j) = HommelGlobal(pselected, tot, 1.0);
                }
                // # nocov end
            }


          }

          if (krest(i,j) == 0) pcomp(i,j) = 1.0;
        }
        
        // Compute family weights
        c(i, 0) = 1;
        for(j = 1; j < nfam; j++) {
          if (method == 1) c(i,j) = c(i,j-1) * (1-ErrorFrac(korig(i,j-1), nperfam, gamma[j-1]));
          if (method == 2) c(i,j) = c(i,j-1) * (1-ErrorFrac(krest(i,j-1), nrest(i,j-1), gamma[j-1]));  
          if (method == 3) c(i,j) = c(i,j-1) * (1-ErrorFrac(krest(i,j-1), nrest(i,j-1), gamma[j-1]));  
        }
        
        // Compute the intersection p-value for the current intersection hypothesis
        pint[i] = IntPvalue(ExtractRow(pcomp, i), ExtractRow(c, i));
        
        // Compute the p-value for each hypothesis within the current intersection
        for(j = 0; j < nhyp; j++) p(i, j) = int_orig(i, j) * pint[i];       
        
      }

    // Compute adjusted p-values
    for (i = 0; i < nhyp; i++) {
        adjpvalue[i] = 0.0;
        temp_vec = ExtractColumn(p, i);
        adjpvalue[i] = *std::max_element(temp_vec.begin(), temp_vec.end());
    }   

    return adjpvalue;

}

vector<double> TradMultAdj(const int &test, const vector<double> &pvalue, const vector<double> &weight, const vector<double> &transition) {

    vector<double> adj_pvalue(pvalue);   

    if (test == 1) adj_pvalue = BonferroniAdj(pvalue, weight);
    // Holm, Hochberg, Hommel
    if (test == 2 || test == 3 || test == 4) adj_pvalue = ClosedTestingAdj(test, pvalue, weight);
    if (test == 5) adj_pvalue = FixedSeqAdj(pvalue, weight);
    if (test == 6) adj_pvalue = ChainAdj(pvalue, weight, transition);

    return adj_pvalue;


}

#endif // MEDDESMULT_H
