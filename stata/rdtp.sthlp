{smcl}
{* *! version 1.0.1  23feb2026}{...}
{vieweralsosee "[R] regress" "help regress"}{...}
{viewerjumpto "Syntax" "rdtp##syntax"}{...}
{viewerjumpto "Description" "rdtp##description"}{...}
{viewerjumpto "Algorithm" "rdtp##algorithm"}{...}
{viewerjumpto "Options" "rdtp##options"}{...}
{viewerjumpto "Stored results" "rdtp##results"}{...}
{viewerjumpto "Examples" "rdtp##examples"}{...}
{viewerjumpto "Limitations" "rdtp##limitations"}{...}
{viewerjumpto "References" "rdtp##references"}{...}
{viewerjumpto "Author" "rdtp##author"}{...}

{title:Title}

{phang}
{bf:rdtp} {hline 2} RD tipping-point search for discrete forcing variables


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:rdtp}
{depvar}
{it:forcingvar}
{ifin}{cmd:,}
{opth by(varname)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opth by(varname)}}unit identifier (required); search is performed
    separately for each level{p_end}
{synopt:{opt band:width(#)}}symmetric bandwidth around each candidate cutoff;
    default is {cmd:75}{p_end}
{synopt:{opt search:range(# #)}}restrict candidates to [{it:min}, {it:max}];
    default is all unique values of {it:forcingvar}{p_end}
{synopt:{opt min:obs(#)}}minimum observations required on each side;
    default is {cmd:10}{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(vcetype)}}variance-covariance estimator passed to
    {helpb regress}; e.g., {cmd:vce(robust)},
    {cmd:vce(cluster} {it:varname}{cmd:)}{p_end}

{syntab:Reporting}
{synopt:{opt noi:sily}}display per-unit progress during the search{p_end}
{synopt:{opt l:evel(#)}}confidence level for regression;
    default is {cmd:level(95)}{p_end}

{syntab:Output}
{synopt:{opt sav:ing(filename)}}save per-unit results to a Stata dataset{p_end}
{synopt:{opt replace}}overwrite existing file when using {cmd:saving()}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:rdtp} searches over candidate cutoff values on a discrete forcing
(running) variable to find, for each unit identified by {opt by()}, the
cutoff that maximizes the R-squared from a linear-spline regression
discontinuity (RD) regression.  This is the first step of a "tipping-point"
analysis: identifying where each unit (e.g., school, county, hospital)
appears to impose a threshold on a test score or other discrete assignment
variable.

{pstd}
The command is designed for settings where the true cutoff is unknown and
varies across units.  By sweeping over all unique values of the forcing
variable within an optional {opt searchrange()}, {cmd:rdtp} implements a
data-driven search akin to the threshold-detection method in
{help rdtp##McEachin2020:McEachin, Domina, and Penner (2020)}.

{pstd}
{cmd:rdtp} performs {it:only} the cutoff search.  It does not estimate
second-stage treatment effects, adjust for multiple testing, or perform
bandwidth selection.  Those steps should be conducted separately once the
cutoffs have been identified.  See
{help rdtp##Lee2010:Lee and Lemieux (2010)} for a survey of RD designs.


{marker algorithm}{...}
{title:Algorithm}

{pstd}
For each unit {it:u} in {opt by()}:

{phang2}
1. Obtain the set of unique values of {it:forcingvar} within
   {opt searchrange()} (or all unique values if unspecified).{p_end}

{phang2}
2. For each candidate cutoff {it:c}:{p_end}

{phang3}
a. Construct:{p_end}
{p 16 18 2}
{it:cut}      = 1({it:forcingvar} {c >=} {it:c}){break}
{it:force}    = {it:forcingvar} {c -} {it:c}{break}
{it:interact} = {it:cut} {c *} {it:force}{p_end}

{phang3}
b. Count observations on each side within {opt bandwidth()}.
   Skip if either side has fewer than {opt minobs()} observations.{p_end}

{phang3}
c. Estimate:{p_end}
{p 16 18 2}
{it:depvar} = {it:{c beta}}{sub:0} + {it:{c beta}}{sub:1}{c ·}{it:cut}
             + {it:{c beta}}{sub:2}{c ·}{it:force}
             + {it:{c beta}}{sub:3}{c ·}{it:interact} + {it:{c epsilon}}{p_end}
{p 16 18 2}
using observations with |{it:force}| {c <=} {opt bandwidth()}.{p_end}

{phang3}
d. Skip if the regression fails, the F-statistic is missing, or
   {it:{c beta}}{sub:1} = 0.{p_end}

{phang3}
e. If R{c 178} > best R{c 178} so far, update the best cutoff for this
   unit.{p_end}

{phang2}
3. Store the results for the best cutoff (or missing values if no valid
   cutoff was found).{p_end}

{pstd}
The R-squared criterion is a heuristic for detecting the cutoff that
produces the sharpest discontinuity.  It does not constitute a formal test.
See {help rdtp##Hansen2000:Hansen (2000)} for the econometric foundations
of threshold detection.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt by(varname)} is required.  It specifies the variable that identifies
units (e.g., schools, districts).  The search is performed separately for
each unique value.  {it:varname} may be numeric or string.
{p_end}

{phang}
{opt bandwidth(#)} sets the symmetric bandwidth around each candidate
cutoff.  Only observations with |{it:forcingvar} {c -} {it:c}| {c <=} {it:#}
are included in the regression.  Default is {cmd:75}.
{p_end}

{phang}
{opt searchrange(# #)} restricts the set of candidate cutoffs to values
of {it:forcingvar} in [{it:min}, {it:max}].  {it:min} < {it:max} is
required.  Endpoints may be integers or non-integers.  If not specified,
all unique values of {it:forcingvar} (within the unit's data) are considered.
{p_end}

{phang}
{opt minobs(#)} sets the minimum number of observations required on each
side of a candidate cutoff (within the bandwidth) for the regression to be
attempted.  Default is {cmd:10}.
{p_end}

{dlgtab:SE/Robust}

{phang}
{opt vce(vcetype)} specifies the variance estimator passed to
{helpb regress}.  Any {it:vcetype} accepted by {cmd:regress} is valid;
e.g., {cmd:vce(robust)}, {cmd:vce(cluster} {it:varname}{cmd:)}.  If not
specified, OLS standard errors are used.
{p_end}

{dlgtab:Reporting}

{phang}
{opt noisily} displays per-unit progress output during the search,
including the selected cutoff, R-squared, and discontinuity estimate
for each unit.
{p_end}

{phang}
{opt level(#)} sets the confidence level for the internal regressions.
Default is {cmd:level(95)} or as set by {helpb set level}.
{p_end}

{dlgtab:Output}

{phang}
{opt saving(filename)} saves the per-unit results to a Stata dataset.
The dataset contains one row per unit with the following variables:
{cmd:unit} (str244), {cmd:cutoff}, {cmd:r2}, {cmd:beta}, {cmd:se},
{cmd:tstat}, {cmd:n_left}, {cmd:n_right}, {cmd:n_total},
{cmd:pred_left}, {cmd:pred_right}.
Units for which no valid cutoff was found have missing values.
{p_end}

{phang}
{opt replace} permits {cmd:saving()} to overwrite an existing file.
{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:rdtp} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n_units)}}number of units searched{p_end}
{synopt:{cmd:r(n_found)}}number of units with a valid cutoff{p_end}
{synopt:{cmd:r(n_skipped)}}number of units with no valid cutoff{p_end}
{synopt:{cmd:r(bandwidth)}}bandwidth used{p_end}
{synopt:{cmd:r(minobs)}}minimum observations per side{p_end}
{synopt:{cmd:r(mean_r2)}}mean R{c 178} across found units{p_end}
{synopt:{cmd:r(mean_beta)}}mean discontinuity estimate{p_end}
{synopt:{cmd:r(mean_tstat)}}mean t-statistic{p_end}
{synopt:{cmd:r(level)}}confidence level used{p_end}
{p2colreset}{...}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}}{cmd:rdtp}{p_end}
{synopt:{cmd:r(cmdline)}}full command as typed{p_end}
{synopt:{cmd:r(depvar)}}dependent variable name{p_end}
{synopt:{cmd:r(forcing)}}forcing variable name{p_end}
{synopt:{cmd:r(byvar)}}by-variable name{p_end}
{synopt:{cmd:r(vce)}}vce specification, if any{p_end}
{synopt:{cmd:r(searchrange)}}search range, if specified{p_end}
{synopt:{cmd:r(saving)}}path of saved results file, if specified{p_end}
{p2colreset}{...}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(results)}}({it:n_units} {c *} 10) matrix of per-unit
    results{p_end}
{p2colreset}{...}

{pstd}
Columns of {cmd:r(results)}:
{cmd:cutoff}, {cmd:r2}, {cmd:beta}, {cmd:se}, {cmd:tstat},
{cmd:n_left}, {cmd:n_right}, {cmd:n_total}, {cmd:pred_left},
{cmd:pred_right}.


{marker examples}{...}
{title:Examples}

    {hline}
{pstd}Self-contained example with synthetic data{p_end}

{pstd}
This example generates a dataset of 5 schools (100 students each), where
each school has a different true cutoff on a discrete test-score variable.
Run this directly from {cmd:help rdtp} to see {cmd:rdtp} recover the
cutoffs.{p_end}

{phang2}{cmd:. * Generate synthetic data for demonstration}{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. set obs 500}{p_end}
{phang2}{cmd:. gen school = ceil(_n / 100)}{p_end}
{phang2}{cmd:. gen test_score = 200 + int(runiform() * 301)}{p_end}
{phang2}{cmd:. * Each school has a different placement cutoff}{p_end}
{phang2}{cmd:. gen true_cutoff = cond(school==1, 300, cond(school==2, 350, cond(school==3, 400, cond(school==4, 325, 375))))}{p_end}
{phang2}{cmd:. gen above = (test_score >= true_cutoff)}{p_end}
{phang2}{cmd:. gen outcome = 50 + 0.1 * test_score + 10 * above + rnormal(0, 5)}{p_end}
{phang2}{cmd:. * Run rdtp to recover the cutoffs}{p_end}
{phang2}{cmd:. rdtp outcome test_score, by(school) searchrange(250 450) noisily}{p_end}

    {hline}
{pstd}Basic usage{p_end}

{phang2}{cmd:. rdtp math_score test_score, by(school_id)}{p_end}

    {hline}
{pstd}Restrict the search to test scores between 225 and 500{p_end}

{phang2}{cmd:. rdtp math_score test_score, by(school_id) searchrange(225 500)}{p_end}

    {hline}
{pstd}Cluster-robust standard errors, narrow bandwidth{p_end}

{phang2}{cmd:. rdtp math_score test_score, by(school_id) bandwidth(50) vce(cluster test_score)}{p_end}

    {hline}
{pstd}Save results and display per-unit progress{p_end}

{phang2}{cmd:. rdtp math_score test_score, by(school_id) saving(rd_results) replace noisily}{p_end}

    {hline}
{pstd}Post-processing the saved results{p_end}

{phang2}{cmd:. use rd_results, clear}{p_end}
{phang2}{cmd:. * Distribution of detected cutoffs}{p_end}
{phang2}{cmd:. histogram cutoff if !missing(cutoff), width(5) title("Detected Cutoffs")}{p_end}
{phang2}{cmd:. * Summary of fit and discontinuity estimates}{p_end}
{phang2}{cmd:. summarize r2 beta tstat if !missing(cutoff), detail}{p_end}

    {hline}
{pstd}Access stored results programmatically{p_end}

{phang2}{cmd:. rdtp math_score test_score, by(school_id) searchrange(225 500)}{p_end}
{phang2}{cmd:. display "Found cutoffs for " r(n_found) " of " r(n_units) " units"}{p_end}
{phang2}{cmd:. matrix list r(results)}{p_end}

    {hline}


{marker limitations}{...}
{title:Limitations}

{phang}
1. {bf:Discrete forcing variable.}  {cmd:rdtp} is designed for settings
where the forcing variable takes on a finite number of distinct values
(e.g., test scores reported as integers).  With continuous forcing
variables, the set of candidates may be extremely large and the search
computationally expensive.
{p_end}

{phang}
2. {bf:R{c 178} is a heuristic.}  Selecting the cutoff that maximizes
R{c 178} does not constitute a formal statistical test for the existence
of a threshold.  The R{c 178} criterion can be influenced by the
distribution of the forcing variable and by heteroskedasticity.
{p_end}

{phang}
3. {bf:No second-stage estimation.}  {cmd:rdtp} identifies candidate
cutoffs but does not estimate causal treatment effects.  After identifying
cutoffs, the researcher should conduct a full RD analysis using appropriate
methods (e.g., {helpb rdrobust} if installed).
{p_end}

{phang}
4. {bf:Symmetric bandwidth.}  The bandwidth is applied symmetrically
around each candidate cutoff.  Asymmetric bandwidths are not supported.
{p_end}

{phang}
5. {bf:No multiple-testing correction.}  When searching over many
candidate cutoffs, the "best" cutoff may be selected by chance.  {cmd:rdtp}
does not adjust for this multiplicity.  Researchers should assess the
robustness of the selected cutoff using domain knowledge and sensitivity
analyses.
{p_end}


{marker references}{...}
{title:References}

{marker Hansen2000}{...}
{phang}
Hansen, B. E.  2000.
Sample splitting and threshold estimation.
{it:Econometrica} 68(3): 575-603.
{browse "https://doi.org/10.1111/1468-0262.00124"}.

{marker Lee2010}{...}
{phang}
Lee, D. S., and T. Lemieux.  2010.
Regression discontinuity designs in economics.
{it:Journal of Economic Literature} 48(2): 281-355.
{browse "https://doi.org/10.1257/jel.48.2.281"}.

{marker McEachin2020}{...}
{phang}
McEachin, A., T. Domina, and A. Penner.  2020.
Heterogeneous effects of early algebra across California middle schools.
{it:Journal of Policy Analysis and Management} 39(3): 772-800.
{browse "https://doi.org/10.1002/pam.22202"}.


{marker author}{...}
{title:Author}

{pstd}
Andrew McEachin


{title:Also see}

{psee}
Online: {helpb regress}, {helpb rdrobust} (if installed)
{p_end}
