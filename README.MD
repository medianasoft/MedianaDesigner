# MedianaDesigner 

<!-- badges: start -->
[![R-CMD-check](https://github.com/medianasoft/MedianaDesignerPrivate/workflows/R-CMD-check/badge.svg)](https://github.com/medianasoft/MedianaDesignerPrivate/actions)
[![codecov](https://codecov.io/gh/medianasoft/MedianaDesignerPrivate/branch/master/graph/badge.svg?token=NwvPou45uf)](https://codecov.io/gh/medianasoft/MedianaDesignerPrivate)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/MedianaDesigner)](https://cran.r-project.org/package=MedianaDesigner)
[![CRAN\_Logs\_Badge](http://cranlogs.r-pkg.org/badges/MedianaDesigner)](https://cran.r-project.org/package=MedianaDesigner)
[![CRAN\_Logs\_Badge\_Total](http://cranlogs.r-pkg.org/badges/grand-total/MedianaDesigner)](https://cran.r-project.org/package=MedianaDesigner)
<!-- badges: end -->

## Efficient Simulation-Based Power and Sample Size Calculations for a Broad Class of Late-Stage Clinical Trials

This package supports efficient simulation-based power and sample size calculations for a broad class of late-stage clinical trials, including Phase II trials, seamless Phase II/III trials and Phase III trials:

* Adaptive trials with data-driven sample size or event count re-estimation (ADSSMod function).
* Adaptive trials with data-driven treatment selection (ADTreatSel function).
* Adaptive trials with data-driven population selection (ADPopSel function).
* Optimal selection of a futility stopping rule (FutRule function).
* Blinded event prediction in event-driven trials (EventPred function).
* Adaptive trials with response-adaptive randomization (ADRand function) [experimental module].
* Traditional trials with multiple objectives (MultAdj function) [experimental module].
* Cluster-randomized trials (ClustRand function) [experimental module].

For more information on this package, visit [Mediana's web site](http://www.mediana.us/free-software). The latest stable version of the package can be downloaded from [CRAN](https://cran.r-project.org/web/packages/MedianaDesigner/index.html). 

The technical manuals with a detailed description of the statistical methodology implemented in each module can be found on Mediana's web site: 

* [Adaptive trials with data-driven sample size or event count re-estimation](http://www.mediana.us/MedianaDesigner/ADSSMod.pdf).
* [Adaptive trials with data-driven treatment selection](http://www.mediana.us/MedianaDesigner/ADTreatSel.pdf).
* [Adaptive trials with data-driven population selection](http://www.mediana.us/MedianaDesigner/ADPopSel.pdf).
* [Optimal selection of a futility stopping rule](http://www.mediana.us/MedianaDesigner/FutRule.pdf).
* [Blinded event prediction in event-driven trials](http://www.mediana.us/MedianaDesigner/EventPred.pdf).
* [Adaptive trials with response-adaptive randomization](http://www.mediana.us/MedianaDesigner/ADRand.pdf).
* [Traditional trials with multiple objectives](http://www.mediana.us/MedianaDesigner/MultAdj.pdf).
* [Traditional trials with cluster-randomized designs](http://www.mediana.us/MedianaDesigner/ClustRand.pdf).

Additional information and multiple case studies can be found in the [online manual (English)](https://medianasoft.github.io/MedianaDesigner) or [online manual (French)](https://medianasoft.github.io/MedianaDesignerFrench).

An online training course on adaptive designs and clinical trial simulation is available on [Mediana's YouTube channel](https://www.youtube.com/channel/UCAQUAsQxxnzje4rEdJ4870Q).

## Functions

Key functions included in the package:

* ADSSMod: Simulation-based design of adaptive trials with data-driven sample size or event count re-estimation.
* ADSSModApp: Graphical user interface for designing adaptive trials with data-driven sample size or event count re-estimation.
* ADTreatSel: Simulation-based design of adaptive trials with data-driven treatment selection.
* ADTreatSelApp: Graphical user interface for designing adaptive trials with data-driven treatment selection.
* ADPopSel: Simulation-based design of adaptive trials with data-driven population selection.
* ADPopSelApp: Graphical user interface for designing adaptive trials with data-driven population selection.
* FutRule: Simulation-based selection of an optimal futility stopping rule at an interim analysis.
* FutRuleApp: Graphical user interface for an optimal selection of a futility stopping rule.
* EventPred: Simulation-based event prediction in trials with an event-driven design.
* EventPredApp: Graphical user interface for event prediction.
* ADRand: Simulation-based design of adaptive trials with response-adaptive randomization.
* ADRandApp: Graphical user interface for designing adaptive trials with response-adaptive randomization.
* MultAdj: Simulation-based design of traditional trials with multiple objectives.
* MultAdjApp: Graphical user interface for power calculations in traditional trials with multiple objectives.
* ClustRand: Simulation-based design of cluster-randomized trials.
* ClustRandApp: Graphical user interface for power calculations in cluster-randomized trials.
* GenerateReport: Simulation report for any of the modules.

## Feedback

Please let us know if you have questions about the R package or underlying methodology. You can contact us at info at mediana.us.  In addition, you could submit questions, issues or feature requests on the [issues page](https://github.com/medianasoft/MedianaDesigner/issues).

We would like to thank multiple individuals, including Thomas Brechenmacher (IQVIA), Douglas McNair (Gates Foundation) and Thomas Peppard (Certara), for their feedback that helped us improve the package and add new features.

## References

* Ahn, C., Heo, M., Zhang, S. (2015). Sample Size Calculations for Clustered and Longitudinal Outcomes in Clinical Research. Chapman and Hall/CRC.
* Alosh, M., Bretz, F., Huque, M. (2014). Advanced multiplicity adjustment methods in clinical trials. Statistics in Medicine. 33, 693-713.
* Bagiella, E., Heitjan, D.F. (2001). Predicting analysis times in randomized clinical trials. Statistics in Medicine. 20, 2055-2063.
* Chuang-Stein, C., Kirby, S., French, J., Kowalski, K., Marshall, S., Smith, M. K. (2011). A quantitative approach for making go/no-go decisions in drug development. Drug Information Journal. 45, 187-202.
* Dmitrienko, A., Bretz, F., Westfall, P.H., et al. (2009). Multiple testing methodology.  Multiple testing problems in pharmaceutical statistics. Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors). New York: Chapman and Hall/CRC Press.
* Dmitrienko, A., Tamhane, A.C. (2011). Mixtures of multiple testing procedures for gatekeeping applications in clinical trials. Statistics in Medicine. 30, 1473-1488.
* Dmitrienko, A., D'Agostino, R. Sr. (2013). Traditional multiplicity adjustment methods in clinical trials. Statistics in Medicine. 32, 5172-5218.
* Dmitrienko, A., Kordzakhia, G., Brechenmacher, T. (2016). Mixture-based gatekeeping procedures for multiplicity problems with multiple sequences of hypotheses. Journal of Biopharmaceutical Statistics. 26, 758-780.
* Dmitrienko, A., Paux, G. (2017). Subgroup analysis in clinical trials. Clinical Trial Optimization Using R. Dmitrienko, A., Pulkstenis, E. (editors). Chapman and Hall/CRC Press, New York.
* Dmitrienko, A., D'Agostino, R.B. (2018). Multiplicity considerations in clinical trials. New England Journal of Medicine. 378, 2115-2122.
* Hayes, R.J., Moulton, L.H. (2009). Cluster Randomised Trials: A Practical Approach. Chapman and Hall/CRC.
* Herson, J. (1979). Predictive probability early termination plans for Phase II clinical trials. Biometrics. 35, 775-783.
* Kordzakhia, G., Brechenmacher, T., Ishida, E., Dmitrienko, A., Zheng, W.W., Lie, D.F. (2018). An enhanced mixture method for constructing gatekeeping procedures in clinical trials. Journal of Biopharmaceutical Statistics. 28, 113-128.
* Millen, B., Dmitrienko, A., Ruberg, S., Shen, L. (2012). A statistical framework for decision making in confirmatory multi-population tailoring clinical trials. Drug Information Journal. 46, 647-656.
* Wang, D., Cui, L., Zhang, L., Yang, B. (2014). An ROC approach to evaluate interim go/no-go decision-making quality with application to futility stopping in the clinical trial designs. New Developments in Statistical Modeling, Inference and Application. Jin, Z., Liu, M., Luo, X. (editors). Springer, New York. 121-147.
* Wassmer, G., Brannath, W. (2016). Group Sequential and Confirmatory Adaptive Designs in Clinical Trials. New York: Springer.
