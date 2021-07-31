# MedianaDesigner 

<!-- badges: start -->
[![R-CMD-check](https://github.com/medianasoft/MedianaDesignerPrivate/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/medianasoft/MedianaDesignerPrivate/actions)
[![codecov](https://codecov.io/gh/medianasoft/MedianaDesignerPrivate/branch/master/graph/badge.svg?token=NwvPou45uf)](https://codecov.io/gh/medianasoft/MedianaDesignerPrivate)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/MedianaDesigner)](https://cran.r-project.org/package=MedianaDesigner)
[![CRAN\_Logs\_Badge](http://cranlogs.r-pkg.org/badges/MedianaDesigner)](https://cran.r-project.org/package=MedianaDesigner)
[![CRAN\_Logs\_Badge\_Total](http://cranlogs.r-pkg.org/badges/grand-total/MedianaDesigner)](https://cran.r-project.org/package=MedianaDesigner)
<!-- badges: end -->

## Efficient Simulation-Based Power and Sample Size Calculations for a Broad Class of Late-Stage Clinical Trials

This package provides efficient simulation-based tools for adaptive designs commonly used in Phase III and seamless Phase II/III trials. The following modules are currently included in the package:

* Adaptive designs with data-driven sample size or event count re-estimation (ADSSMod function).
* Adaptive designs with data-driven treatment selection (ADTreatSel function).
* Adaptive designs with data-driven population selection (ADPopSel function).
* Optimal selection of a futility stopping rule (FutRule function).
* Event prediction in event-driven trials (EventPred function).

For more information on this package, visit [Mediana's web site](http://www.mediana.us/free-software). The latest stable version of the package can be downloaded from [CRAN](https://cran.r-project.org/web/packages/MedianaDesigner/index.html). 

The technical manuals with a detailed description of the statistical methodology implemented in each module can also be found on Mediana's web site: 

* [Adaptive designs with data-driven sample size or event count re-estimation](http://www.mediana.us/MedianaDesigner/ADSSMod.pdf).
* [Adaptive designs with data-driven treatment selection](http://www.mediana.us/MedianaDesigner/ADTreatSel.pdf).
* [Adaptive designs with data-driven population selection](http://www.mediana.us/MedianaDesigner/ADPopSel.pdf).
* [Optimal selection of a futility stopping rule](http://www.mediana.us/MedianaDesigner/FutRule.pdf).
* [Event prediction in event-driven trials](http://www.mediana.us/MedianaDesigner/EventPred.pdf).

Additional information and multiple case studies can be found in the [online manual](https://medianasoft.github.io/MedianaDesigner).

## Functions

Key functions included in the package:

* ADSSMod: Simulation-based design of adaptive trials with data-driven sample size or event count re-estimation.
* ADSSModApp: Graphical user interface to design adaptive trials with data-driven sample size or event count re-estimation.
* ADTreatSel: Simulation-based design of adaptive trials with data-driven treatment selection.
* ADTreatSelApp: Graphical user interface to design adaptive trials with data-driven treatment selection.
* ADPopSel: Simulation-based design of adaptive trials with data-driven population selection.
* ADPopSelApp: Graphical user interface to design adaptive trials with data-driven population selection.
* FutRule: Simulation-based selection of an optimal futility stopping rule at an interim analysis.
* FutRuleApp: Graphical user interface for an optimal selection of a futility stopping rule.
* EventPred: Simulation-based event prediction in trials with an event-driven design.
* EventPredApp: Graphical user interface for event prediction.
* GenerateReport: Simulation report for any of the modules.

## References

* Bagiella, E., Heitjan, D.F. (2001). Predicting analysis times in randomized clinical trials. Statistics in Medicine. 20, 2055-2063.
* Chuang-Stein, C., Kirby, S., French, J., Kowalski, K., Marshall, S., Smith, M. K. (2011). A quantitative approach for making go/no-go decisions in drug development. Drug Information Journal. 45, 187-202.
* Dmitrienko, A., Paux, G. (2017). Subgroup analysis in clinical trials. Clinical Trial Optimization Using R. Dmitrienko, A., Pulkstenis, E. (editors). Chapman and Hall/CRC Press, New York.
* Dmitrienko, A., D'Agostino, R.B. (2018). Multiplicity considerations in clinical trials. New England Journal of Medicine. 378, 2115-2122.
* Herson, J. (1979). Predictive probability early termination plans for Phase II clinical trials. Biometrics. 35, 775-783.
* Millen, B., Dmitrienko, A., Ruberg, S., Shen, L. (2012). A statistical framework for decision making in confirmatory multi-population tailoring clinical trials. Drug Information Journal. 46, 647-656.
* Wang, D., Cui, L., Zhang, L., Yang, B. (2014). An ROC approach to evaluate interim go/no-go decision-making quality with application to futility stopping in the clinical trial designs. New Developments in Statistical Modeling, Inference and Application. Jin, Z., Liu, M., Luo, X. (editors). Springer, New York. 121-147.
* Wassmer, G., Brannath, W. (2016). Group Sequential and Confirmatory Adaptive Designs in Clinical Trials. New York: Springer.