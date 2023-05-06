This code replicates the empirical exercise in Bertsche, Braun (2020):
Identification of Structural Vector Autoregressions by Stochastic Volatility

If you run "runme.m", the code starts to replicate the results from the empirical exercise by:
1) calling "Main_Estimation.m", which estimates the oil market SVAR with 12 lags and 4 heteroskedastic shocks.
% The results are saved in the output folder ("results_bh19.mat"). This may take a bit of time.
2) calling "Main_Testing_Identification.m", which will estimate the SVAR model for 
% 0,1,2 and 3 heteroskedastic shocks in order to test for the heteroskedasticity rank.  This may take a bit of time.
% Results are saved in output folder ("results_bh19_identification.mat")
3) calling "Application_Table_1_Information_Criteria.m", which prints Table 1. 
% Note that we do not provide code for the other models. 
4) calling "Application_Table_2_TestOfIdentification.m", which prints Table 2.
5) calling "Application_Table_3_TestingIVs.m", which computes the overidentification tests and prints Table 3
5) calling "Application_Figure_1_IRFs.m", which computes Figure 1 (Impulse Response Functions.
% Note that for the SVAR-IV we use code from  Montiel Olea, Stock and Watson (2020, JOE forthcoming) 
6) calling "Application_Table_4_FEVD_RPO.m", which prints Table 4.

We also prepared a smaller example, which relies on simulated data ("Simple_Example.m")
This code simulates data from a bivariate SV-SVAR and compares estimates and confidence intervals from EM-1 and EM-2.
As you will see: at least for this simple model, EM-1 will be much faster and the results are virtually the same.