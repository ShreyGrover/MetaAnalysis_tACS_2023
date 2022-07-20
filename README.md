# Computation of effect sizes 

### Compute_ES.R

Computes effect sizes using raw data

**1. Computing effect sizes using raw data (means, SDs, t-statistics)**

functions:

`compute_G_within()`  
Computes Hedges’ G in within-subjects experiments

	n: sample size
	M1: active mean
	M2: sham mean
	SD1: active standard deviation
	SD2: sham standard deviation
	r: correlation between conditions 

`compute_G_between()`  
Computes Hedges’ G in between-subjects experiments

	n1: active group size
	n2: sham group size
	M1: active mean
	M2: sham mean
	SD1: active standard deviation
	SD2: sham standard deviation

`compute_G_diffScore()`  
Computes Hedges’ G in within-subjects experiments using difference score and its SD

	n: sample size
	Mdiff: difference score (difference between active and sham)
	SDdiff: standard deviation of the difference score
	r: correlation between conditions

`compute_G_indT()`  
Computes Hedges’ G in between-subjects experiments using independent samples t-test statistic
	
	n1: active group size
	n2: sham group size
	t: independent samples t-test statistic

**2. Combining effects**

The `comb_ef` column in the `metadata` data frame has markers for effects to be combined. The script uses Borenstein et al. book formulas to compute average effects and their respective SDs. 

# Pooling effect sizes 

## Outcome-based analysis

### Meta_Effects_Outcome_Analysis.RmD

Pools effect sizes for All, Performance and RT
Computes effects for cognitive domains
Older adults and clinical subgroups analysis
Publication Bias

**1. Pooling effect sizes**

`robu()` function: 
- for pooled effect sizes with N<40 use small = TRUE, which turns on small sample adjustments (see Tipton, 2015 for details)
- for regression always recommended to use small = TRUE (Tipton, 2015)

`sensitivity()` function takes in the model output from the `robu()` function and runs sensitivity analysis by varying correlation between the effects nested within a study using r = {0, 0.2, 0.4, 0.6, 0.8, 1}. 

Robumeta documentation: https://cran.r-project.org/web/packages/robumeta/robumeta.pdf

**2. Outliers removal**

Experiments with absolute values of studentized deleted residuals above 1.96 are considered outliers (Viechtbauer & Cheung, 2010)

`rstudent()` function (metafor package) computes studentized deleted residuals
`rstudent()` accepts only metafor-produced objects, so we are using metafor function `rma.uni()` which performs univariate random effects meta-analysis.  Important note: we use it with the effect weights derived from fitting RVE model first. 

Outlier removal function:

 `remove_rst()`

	df: data set (in case of our meta-analysis “All”, “Performance”, or “RT”	
	model_domain: Domain-specific RVE model fitted using robu() function
	domain: name of the cognitive domain (string)
   

**3. Forest plotting function** 
(not the one we use for figures, but useful for quickly displaying the data)

`forest_rve()`

 	model_intercept: RVE model fitted using robu() function
	domain: name of the cognitive domain (string)
	scale_size: scaling factor for effect size markers. The markers represent effect size 	weights computed by robu(). For display purposes they are additionally scaled by 	scale_size, as the markers might be too small or too big deducting readability of the plot.

**4. Publication bias**

Uses several metafor functions:

* `rma.uni()` - univariate random effects meta-analysis, we use it with the effect weights derived from fitting RVE model first
* `funnel()` - funnel plotting function
* `regtest.rma()` - Egger’s regression (we are not using it anymore, kept in script in case needed)
* `trimfill()` - runs trim and fill procedure and returns adjusted effect

Detailed documentation can be looked up here: https://cran.r-project.org/web/packages/metafor/metafor.pdf

## Hypothesis-based analysis

### Meta_Effects_Hypothesis_Analysis.RmD

Pools effect sizes for All, Performance and RT
Computes effects for cognitive domains (no outlier removal, but easy to add if needed)
Older adults (no outlier removal)

## Phase subgroup analysis 
### Meta_Subgroup_Phase.Rmd
Pools effect sizes for All, Performance and RT in a subgroup of studies manipulating phase intentionally with a defined hypothesis

**1. Overall effect of phase manipulation**  
Pools effects of tACS in in-phase and anti-phase studies combined.

**2. Effect of in-phase tACS on improvement**  
Pools effects of tACS in a subset of experiments that performed intentional in-phase manipulation in order to improve cognitive or clinical outcome.

**3. Effect of anti-phase tACS on impairment**  
Pools effects of tACS in a subset of experiments that performed intentional anti-phase manipulation in order to disrupt cognitive function. 

# Regressions

### Meta_Regression_Revised.R

**1. Simple Regressions for All, Performance and RT:**

Predictors:

	individual frequency (IF)
	intensity (numeric)
	intensity (categorical/ low (<=1 mA), high (>1 mA))
	duration (categorical/ <= 20 min, > 20 min)
	timing of assessment (categorical/ online, offline)
	timing of stimulation (categorical/ during rest or task)
	HD (categorical/ yes/no)
	current modeling (categorical/ yes/no)
	neuroguided (categorical/ yes/no)
	design (categorical/ within, between-subjects)
	age (categorical/ young, older adults)
	blinding (categorical/ single, double)

**2. Interactions with timing of assessment:**

*2.1 Probing Interactions:*  

Likelihood ratio tests using R function `anova()` (https://www.rdocumentation.org/packages/car/versions/3.0-11/topics/Anova) are used to compare 2 models: with and without the interaction. `anova()` doesn’t accept robumeta objects, so we use metafor function `rma()` with the weights of the effects derived from RVE model fitted using robu() function. 

Probing multiple interactions requires correction for multiplicity, so we divide alpha by the number of interactions we are testing. 

*2.2 Pairwise Comparisons:*  

After confirming that interaction is significant, we perform pairwise comparisons. Here, again, correction for multiplicity is required and is done by dividing alpha by the number of pairwise comparisons.  

To perform comparisons we use `anova()` function with the `rma()` object (again, with weights derived from RVE model). Details of this process are described here: https://www.metafor-project.org/doku.php/tips:multiple_factors_interactions

### Meta_Regression_WM_Revised.Rmd
Simple meta-regressions for All and Performance outcomes in Working Memory domain. 

### Meta_Regression_LTM_Revised_LTM.Rmd
Simple meta-regressions for All and Performance outcomes in Long-Term Memory domain. 
	
### Meta_Regression_OtherDomains_Revised.Rmd
Simple meta-regressions for All and Performance outcomes in Executive Control, Attention, Intelligence, Motor Learning and Motor Memory domains. 
	
# Figures

### Meta_Plots_Fig1_studies_by-years.Rmd
Fig. 1 (saved and HTML and svg)

### Meta_Plots_Fig2_Flowchart.Rmd	
Fig. 2 (saved and HTML and svg)

### Meta_Plots_Fig4_Treemap.Rmd
Fig. 4 (saved and HTML and svg)

### Meta_Plots_FigS1_Treemap.Rmd
Fig. S1 (saved and HTML and svg)

### Meta_Plots_Fig5-8_S2.Rmd
Figs. 5,6,7,8,S2 (saved and HTML and svg)

### Meta_Plots_FigS4_ROB.Rmd
Supplementary Fig. 4 (saved and HTML and svg)


# Data Files

### Raw data
**(1) ES_raw_hypothesis.xlsx**  

**(2) ES_raw_outcome.xlsx**  

**(3) ES_raw_NEW_EFFECTS_nov21.xlsx**  

These files contain means and standard deviations as well all other data collected form the included studies. (1) contains studies that will be included in the hypothesis-based analysis, (2) contains studies that will be included in the outcome-based analysis, and (3) contains new studies added in November 2021 after additional search through the bibliographies (*Note: (3) includes studies that will be included in either hypothesis- or outcome-based analysis*). These files are used in the script `Compute_ES.R`, which takes in the raw data and computes Hedges’ G and its variance for each experiment. Column `comb_ef` is used for averaging effects: same numerical value is entered for the effects to be averaged (same for a group of effects, but different for different groups of effects), and 0 is entered otherwise. Output of Compute_ES.R is saved in intermediate data files `ES_Computed_Outcome.xlsx`, `ES_Computed_Hypothesis.xlsx` and `ES_Computed_NEW.xlsx` which are further used to create/update the files used for analysis (see below). 

### Data for analysis 

**Analysis_Outcome_Revised.xlsx**  

Contains all data for outcome-based analysis, i.e. effects from experiments aiming to improve the outcome or exploratory in nature. The data in this file was constructed using `ES_Computed_Outcome.xlsx` and outcome-based effects from `ES_Computed_NEW.xlsx` (outputs of `Compute_ES.R`).  

*Note on effect size columns:*  

- `G_raw`: Hedges’ G as computed by Compute_ES.R
- `seG`: standard error of Hedges’ G as computed by Compute_ES.R
- `reversal`: values either equal to G_raw or are reversed in sign based on the expected direction of the effects guided by the goal of improvement of the outcome
- `y`: this column is a copy of reversal column (with a formula). This is the column used in all analyses as Hedges’ G. 
- `v`: squared seG. This is the column used in all analyses as variance of G.  
 
***Important: don’t change the names y and v as metafor functions recognize them by default as effect size and its variance.*** 

*Note on the `ID` column:*  

Contains the ID of the study that is the same for all experiments nested within a study. This column values are used by `robu()` to identify groupings of effects within a study and adjust their weights. The dataset is sorted alphabetically by the name of the experiment (column `experiment`), so If new studies are added to the dataset, ID column has to be adjusted: IDs assigned to new studies while maintaining alphabetical order. 


**Analysis_Hypothesis_Revised.xlsx**  

Contains all data for hypothesis-based analysis, i.e. effects from experiments aiming to improve or impair the outcome. The data in this file was constructed using `ES_Computed_Hypothesis.xlsx` and hypothesis-based effects from `ES_Computed_NEW.xlsx` (outputs of `Compute_ES.R`).  

*Note on effect size columns:*  

- `G_raw`: Hedges’ G as computed by Compute_ES.R
- `seG`: standard error of Hedges’ G as computed by Compute_ES.R
- `reversal`: values either equal to G_raw or are reversed in sign based on expected direction of the effects guided by the hypothesis  
- `y`: this column is a copy of reversal column (with a formula). This is the column used in all analyses as Hedges’ G. 
- `v`: squared seG. This is the column used in all analyses as variance of G.  
 
***Important: don’t change the names y and v as metafor functions recognize them by default as effect size and its variance.*** 

**New_effects_Nov2021.xlsx**  

Separate data file that contains new studies added in November 2021. Not used in any analysis, kept for the record of new studies. Use only for quick reference of a new study; all the new studies are already added to the main analysis datasets (`Analysis_Outcome.xlsx` and `Analysis_Hypothesis.xlsx`).

**Treemap_OutcomeBasedEffects.xlsx**  

This dataset is used for generating the treemap figure in `MA_Treemap.Rmd`. Dataset is similar to Appendix 1, except for it has separate columns for author’s name, experiment id, and year for plotting purposes and only consists of effects for exploratory experiments or experiments hypothesizing an improvement in cognitive function.

**Treemap_OutcomeBasedEffects.xlsx**  

This dataset is used for generating the treemap figure in `MA_Treemap.Rmd`. Dataset is similar to Appendix 1, except for it has separate columns for author’s name, experiment id, and year for plotting purposes and only consists of effects for experiments hypothesizing a disruption in cognitive function.

**robdata.xlsx**  

The data file contains RoB assessments and is used to generate the RoB plot in `MA_ROB.Rmd`. 

**Search_results_data.xlsx**  

Results of study search and selection. Not used in R scripts and notebooks. Kept for the record and for manual input of `values` into `MA_Flowchart.Rmd`.

**Studies_by_years.xlsx**  

List of individual studies included in meta-analyses and their years of publication. Used for generating the figure in `MA_studies_by-years.Rmd`.

# Effect Size Reversals  

### Outcome-based analysis  

![](https://github.com/renatafay/tACS_meta_analysis_2021/blob/main/effects_reversal_outcome.png)  
  
### Hypothesis-based analysis  

![](https://github.com/renatafay/tACS_meta_analysis_2021/blob/main/effects_reversal_hypothesis.png)  


