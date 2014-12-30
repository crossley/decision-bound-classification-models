This package of files contains model fitting programs for classification of 2 dimensional stimuli as described by General Recognition Theory (Ashby & Townsend 1986). The core style of this code was written by John Ennis way back circa 2006.  I have modified it heavily since then, and so the current result is a bit of a stylistic melting pot. It nevertheless gets the job done. 

The package contains the following rule types:

Unidimensional X
Unidimensional Y
General Conjunctive Classifier - Unequal Variance
General Conjunctive Classifier - Equal Variance
General Linear Classifier
General Quadratic Classifier
Hybrid - Unequal Variance
Hybrid - Equal Variance
Hybrid - Lin Only
Hybrid - Optimal

The model fitting program fits the models with both possible category labels and returns the best fit of the two.  

The hybrid and conjunctive models expect one category to be in the upper left only.

**** Explanation of folders ****

"nll_prog" contains all of the programs that compute the negative log likelihood of a certain pattern of responses given a specific rule.
	Each of these programs is a function that has the decision bound parameters as its input and the nll as output.  In each case, the variable
	"data" is globally defined.

"data" cotains the raw (preprocessed) data as well as the procesed data.  The raw data can be in any form, but the user must specify in the program
	"transform_data" how to process the data into the correct form (See "transform_data" explanation).  "data" also contains the files of user
	specified starting bounds for the conjunctive and hybrid rules in "fitting_script."

"output" contains the results of "fitting_script."

"servant_progs" contains the various helper scripts called by either "fitting_script" or "data_analysis_gui."  Within "servant_progs" there are several 
	subfolders that contain possibly outdated helper scripts.

**** Explantion of individual programs ****

"fitting_script" performs the macro level analysis whereby many data sets can be analyzed using all of the rules.  "fitting_script" allows the user to 
	specify the subject numbers in an array, to tell the program whether to fit the responses or the true categories, which hybrid models to 
	fit if any, as well as the number of trials, the number of blocks and the length of each block.
	
	"fitting_script" calls "transform_data" to process the data, "record_starting_bounds" to allow the user to specify starting bounds for the gcc
	and hybrid rules, and "run_model_fits" to perform the actual model fitting.
	
"transform_data" MUST be written by the user to transform the data into the form [corr_cat x y resp].  The values of x and y should be numbers between
	0 and 100, and if an optimal hybrid bound is to be fit its unidimensional x bound should be at x = 50.  The categories should be divided so that 
	if a gcc is fit one category is in the upper left (second quadrant).  The output of "transform_data" is expected to be subject###.dat numbered
	by subject.
	
"record_starting_bounds" calls the servant scripts that load the user specified starting bounds for the gcc and hybrid models.

"run_model_fits" runs the actual model fitting.  It produces a master output file called "all_results.dat" that contains a listing of the best
	fitting models for each data set, with fit measured either by BIC or nll.  "run_model_fits" also produces an output file for each subject
	labeled "sub####Block#_results.dat" in which model by model fit information is given for each subject.
	
	



