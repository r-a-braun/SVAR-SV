clear; clc;
load('output/results_bh19.mat')
FEVD_RPO = zeros(specs.hor,2); 
horizons = [0:24:72]; horizons(1) = 12; 
fevds = [ IRF.fevd(horizons,(1*K-K)+3 ), IRF.fevd(horizons,(3*K-K)+3 ) ];
lq = [ IRF.fevd_ql(horizons,(1*K-K)+3 ), IRF.fevd_ql(horizons,(3*K-K)+3 ) ];
uq = [ IRF.fevd_qu(horizons,(1*K-K)+3 ), IRF.fevd_qu(horizons,(3*K-K)+3 ) ];  
Data = 100.*[fevds , lq,uq];
Data = round(Data(:,[1,3,5,2,4,6])',2);

 
%% Print FEVD shares  
input.data = Data;
input.tableColLabels = {' h=1', ' h=24', ' h=48', ' h=72'};
input.tableRowLabels = {'$\varepsilon_{1t}$', '5% quantile', '95% quantile', '$\varepsilon_{3t}$', '5% quantile', '95% quantile'};
input.dataFormat = {'%.2f', 4};  
input.tableCaption = 'Forecast Error Variance Decomposition Real Oil Price';
input.tableLabel = 'FEVDs';
FEVDtable = latexTable(input); 


