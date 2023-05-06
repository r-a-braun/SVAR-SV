clear; clc;
load('output/results_bh19.mat') 
load('data/IVs_oilsupply.mat') 

%% Run the Tests 
Ehat = Param_est.e;  
IVS = [IVs.K09,IVs.BH19,IVs.K08a, IVs.CCI,IVs.DK19]; 
startdate = datetime(1974,1,1);  
idx_Y = datenum(date)>=datenum(startdate); 
IVS = IVS(idx_Y,:); 
ivs_label = {'K09', 'BH19', 'K08', 'CCI19', 'DK19'};  
IVS = IVS(specs.p +1:end,:);  
input.param_est = Param_est;
input.VarTheta = VarTheta;
input.m = IVS;
input.specs = specs;
[ output ] = test_IV_2step( y,input );


%% Print Results: Part 1) Correlations 
input.data = output.corr;
input.tableColLabels = ivs_label;
input.tableRowLabels = {'$\varepsilon_{1t}$','$\varepsilon_{2t}$',...
    '$\varepsilon_{3t}$','$\varepsilon_{4t}$'};
input.dataFormat = {'%.2f',5}; 
input.tableCaption = 'Regression Coefficients';
input.tableLabel = 'TestsIV';
correlations = latexTable(input); 

%% Print Results: Part 2) Regression coefficients with s.e.
input.data = [output.phi;output.se];
input.tableColLabels = ivs_label;
input.tableRowLabels = {'$\psi_{1}$','$\psi_{2}$',...
    '$\psi_{3}$','$\psi_{4}$',...
    '$\psi_{1} (s.e.) $','$\psi_{2} (s.e.) $',...
    '$\psi_{3} (s.e.) $','$\psi_{4} (s.e.) $'};
input.dataFormat = {'%.2f',5};  
input.tableCaption = 'Correlations';
input.tableLabel = 'TestsIV';
coefficients = latexTable(input); 


%% Print Results: Part 3) test results  
input.data = [output.indices; output.nonrelevance; output.exogeneity];
input.tableColLabels = ivs_label;
input.tableRowLabels = {'shock','relevance','exogeneity'};
input.dataFormat = {'%.2f',5};  
input.tableCaption = 'Test for relevance and exogeneity';
input.tableLabel = 'TestsIV';
tests = latexTable(input); 

 
