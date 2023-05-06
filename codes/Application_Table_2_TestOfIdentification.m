
clear; clc;
load('output/results_bh19_identification.mat')
%% Print Test Results (h=1)
input.data = [Q1_stat(:,1), Q1_dof(:,1), Q1_p(:,1),...
    Q2_stat(:,1), Q2_dof(:,1), Q2_p(:,1),...
    LM_stat(:,1), LM_dof(:,1), LM_p(:,1)];
input.tableColLabels = {'$Q_1(h)$', 'df', '$p$-val', '$Q_2(1)$', 'df', '$p$-val', '$LM(1)$', 'df', '$p$-val'};
input.tableRowLabels = {'$r_0=0$','$r_0=1$','$r_0=2$','$r_0=3$'};
input.dataFormat = {'%.2f',9};  
input.tableCaption = 'Test for Identification';
input.tableLabel = 'Teststats';
test_identification = latexTable(input); 