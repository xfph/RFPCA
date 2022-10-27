# This file contains a brief introduction to the code package.

The folder prog contains function files of the proposed RFPCA and other related methods in the paper, detailed as follow.
mvt.m
-- mvt obtains the ML estimate for the Mt distribution.

mvn.m
-- mvn obtains the ML estimate for the matrix-normal distribution.

mxvt.m
-- mxvt obtains the ML estimate for the MxT distribution.

mt.m
-- mt obtains the ML estimate for the multivariate t distribution.

xxx_ini.m
-- xxx_ini initializes the parameters for method xxx.

bpca.m
-- bpca implements bidirectional PCA (BPCA).


The folder expr contains the codes for sections 4.1 and 4.3.2 in the paper, detailed as follow.

The subfolder MTS contains the codes for the classification of the multivariate time series datasets.
AUS_without.mat, AUS_caseIV.mat
-- The training X and test Y sets are obtained by one random splitting of the AUSLAN dataset. 
    'without' means without outliers; 'caseIV' means the case IV mentioned in Section 4.3, that is, the outlying observations from the uniform distribution U(0,10), and its proportion is 20%.

demo_err_aus.m
-- demo_err_aus shows the lowest error rates of different methods over one random splitting of the AUSLAN dataset. 
    Note: To run the results on the ECG dataset, you can download ECG data from http://www.mustafabaydogan.com./files/viewcategory/20-data-sets.html, and replace the AUSLAN data set in this M-file, by dividing it into training and test sets.

mts_1.m
-- mts_1 calculates the classification error rates of different methods.

compr.m
-- compr obtains compressed representations of different methods in the paper.

nnerr.m
-- nnerr calculates the classification error rates using the 1-nearest neighbor classifier.

test_MxT.m
-- test_MxT shows the results of running the MxT with parsimonious covariance structures on the AUSLAN dataset.

The subfolder simulation contains the codes for all the experiments in section 4.1 of the paper, detailed as follow.
expr_4_1_1.m, expr_4_1_2.m, expr_4_1_3.m, expr_4_1_4.m and expr_4_1_5.m
-- These M-files obtain the results of Sec 4.1.1, Sec 4.1.2, Sec 4.1.3, Sec 4.1.4 and Sec 4.1.5, respectively. 

demo_expr_4_1_1.m
--  demo_expr_4_1_1 produces Figure 1 in the main text. 

demo_expr_4_1_2.m
--  demo_expr_4_1_2 produces Figure 2 in the main text. 

demo_expr_4_1_3.m
-- demo_expr_4_1_3 shows Table 2 in the main text. 

demo_expr_4_1_4.m
--  demo_expr_4_1_4 produces Figure 4 in the main text. 

demo_expr_4_1_5.m
--  demo_expr_4_1_5 produces Figure 5 in the main text. 

The subfolder data contains the synthetic data used for some cases in each simulation experiment.

The subfolder result contains the results of each experiment, detailed as follow.
simu1_low.mat, simu1_high.mat
-- Saved results for Section 4.1.1.

simu2_test_llh_1.mat, simu2_test_llh_2.mat
-- Saved results for Section 4.1.2.

relED_50_1_case_x_xx.mat, ..., relED_50_6_case_x_xx.mat
-- Saved results for Section 4.1.3.

simu4_caseoc1.mat, simu4_caseoc2.mat, simu4_caseoc3.mat
-- Saved results for Section 4.1.4.

simu5_t_itnum_oc.mat
-- Saved results for Section 4.1.5.
