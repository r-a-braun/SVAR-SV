clear; clc;
load('data/BH19.mat')
addpath('functions_svsvar')
y = [100*diff(log(dataBH19.production)),...
    100*diff(log(dataBH19.OECD)),...
    100*log(dataBH19.realWTI(2:end)),...
    dataBH19.dInvent(2:end)];
date = dataBH19.date(2:end); 
startdate = datetime(1974,1,1);  
idx_Y = datenum(date)>=datenum(startdate); 
y = y(idx_Y,:); 
[T,K]=size(y);


%% Specifications
specs.param_init = []; 
p = 12;
specs.p = p;  % lag length
specs.max_it = 500;  % Maximum # iterations
specs.M = 5000; % Draws Importance Sampler 
specs.tol = 1e-03; % Tolerance EM algorithm
specs.alpha = [.1,.32]; % Confidence level IRFs
specs.hor = 72; % Horizons IRF
specs.S_b = eye(K^2); % selection matrix B elements (unrestricted)
specs.scale = 20; % offset constant : std(error)./scale
specs.r = K;
specs.unitshocks  = 1;
specs.M2 = 10000;% Draws Oakes Method 

%% TEST FOR IDENTIFICATION:
tol = 1e-3; max_it = 100;
Q1_stat = zeros(K,3); Q1_p = zeros(K,3);
Q2_stat = zeros(K,3); Q2_p = zeros(K,3);
Q1_dof = zeros(K,3); Q2_dof = zeros(K,3);
LM_dof = zeros(K,3); LM_p = zeros(K,3); LM_stat = zeros(K,3);

star_rep = 10;
for r=0:K-1
    specs.r = r;
    if r==0
        [Bhat, Sigmahat, Uhat,llh_lin] = VARmlike(y,p) ;
        B_r = chol(Sigmahat)';
    else
        llh_nrep = zeros(star_rep,1);
        Uhat_nrep = zeros([size(Uhat),star_rep]);
        Bhat_nrep = zeros(K,K,star_rep);
        for nrep = 1:star_rep
            [Param_est_can , llh_can ]  = SVAR_SV_em_approx( y, specs);
            Uhat_nrep(:,:,nrep) = Param_est_can.u;
            Bhat_nrep(:,:,nrep) = Param_est_can.B;
            llh_nrep(nrep) = llh_can.llh; 
        end 
        [~,idx]=max(llh_nrep); 
        Uhat = Uhat_nrep(:,:,idx);
        B_r = Bhat_nrep(:,:,idx); 
    end
    disp(r)
    [ Q11,  Q21, LM1] = svar_sv_tests( Uhat, B_r, r, 1);
    [ Q13,  Q23, LM3] = svar_sv_tests( Uhat, B_r, r, 3);
    [ Q16,  Q26, LM6] = svar_sv_tests( Uhat, B_r, r, 6);
    Q1_stat(r+1,:) = [Q11.teststat;Q13.teststat;Q16.teststat];
    Q1_p(r+1,:) = [Q11.p;Q13.p;Q16.p];
    Q1_dof(r+1,:) = [Q11.dof;Q13.p;Q16.dof];
    LM_p(r+1,:) = [LM1.p;LM3.p;LM6.p];
    LM_stat(r+1,:) = [LM1.teststat;LM3.teststat;LM6.teststat];
    LM_dof(r+1,:) = [LM1.dof;LM3.dof;LM6.dof];

    Q2_stat(r+1,:) = [Q21.teststat;Q23.teststat;Q26.teststat];
    Q2_p(r+1,:) = [Q21.p;Q23.p;Q26.p];
    Q2_dof(r+1,:) = [Q21.dof;Q23.p;Q26.dof];
end
save('output/results_bh19_identification.mat')