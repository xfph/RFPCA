# This file contains a brief introduction to the code package.

The folder prog contains function files of the proposed RFPCA and other related methods in the paper, details as below.
mvt.m
-- Main file used to obtain the ML estimate for the Mt distribution.

mvn.m
-- Main file used to obtain the ML estimate for the matrix-variate normal distribution.

mxvt.m
-- Main file used to obtain the ML estimate for the T distribution.

mt.m
-- Main file used to obtain the ML estimate for the multivariate t distribution.

xxx_ini.m
-- Code to implement the initialization of the corresponding method parameters.

bpca.m
-- Code to implement bidirectional PCA (BPCA).


The folder expr contains the codes for all the experiments in section 4.1 of the paper, details as below. 

expr_4_1_1.m, expr_4_1_2.m, expr_4_1_3.m, expr_4_1_4.m and expr_4_1_5.m
-- Codes used to obtain the results of Sec 4.1.1, Sec 4.1.2, Sec 4.1.3, Sec 4.1.4 and Sec 4.1.5, respectively. 

demo_expr_4_1_1.m
--  Code to produce Figure 1 in the main text. 

demo_expr_4_1_2.m
--  Code to produce Figure 2 in the main text. 

demo_expr_4_1_3.m
--  Code to show Table 2 in the main text. 

demo_expr_4_1_4.m
--  Code to produce Figure 4 in the main text. 

demo_expr_4_1_5.m
--  Code to produce Figure 5 in the main text. 

The subfolder data contains synthetic data for some cases of each simulation experiment.

The subfolder result contains the relevant results of each experiment, with details below. 
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




