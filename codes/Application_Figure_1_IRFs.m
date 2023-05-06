clear; clc;
load('output/results_bh19.mat')
load('data/IVs_oilsupply.mat')  

warning('off')
addpath('functions_svsvar')
addpath(genpath('functions_olea_et_al'))

Ehat = Param_est.e;
IVS = [IVs.K09,IVs.BH19,IVs.K08a, IVs.CCI,IVs.DK19];
startdate = datetime(1974,1,1);
idx_Y = datenum(date)>=datenum(startdate);
IVS = IVS(idx_Y,:);
ivs_label = {'K09', 'BH19','K08','CCI19','DK19'};
dataset_name =  'oilvar' ;
[Bhat, Sigmahat, Uhat,llh_lin] = VARmlike(y,specs.p) ;
idx_IVs = [1,2,5];
idx_Unit = [1,1,3];
direct = pwd;
savdir = strcat(direct,'/figures/');  %selected directory where the output files will be saved
horizon = (0:specs.hor)';
p = specs.p;
confidence  = 1-specs.alpha(1);
columnnames = [{' prod' },  {' Activity'}, {' price'}, {' invent'}];
time        = 'Month';
NWlags      = 0;
scale       = 1;
horizons    = specs.hor;
IRFselect   = [];
cumselect = [];
a = 1;
variables = {' q_t',  ' wip_t', ' p_t', ' i_t'};
shocks=  {' \varepsilon_{1t}', ' \varepsilon_{2t}',' \varepsilon_{3t}',' \varepsilon_{4t}',};
a2 = 1;
csdummy = [1,1,0,1];
hFig1 = figure(1);
negative = [1,1,0];
for a = 1:size(idx_IVs,2)
    z = IVS(:,idx_IVs(a));  z(isnan(z)==1)=0;
    ydata = y; 
    norm   = idx_Unit(a);
    [Plugin, InferenceMSW] = SVARIV(p, confidence, ydata, z, NWlags, norm, scale, horizons, savdir, columnnames, IRFselect, cumselect, time, dataset_name);
    ii = idx_Unit(a); 
    for i = 1:K
        subplot(size(idx_IVs,2),K,a2)
        if csdummy(i)==1
            if negative(a)==1
                hold on
                plot(horizon, -IRF.irfs_cs(:,(ii*K-K)+i ),'k-', ...
                    horizon,-[squeeze(IRF.IRF_cs_qu( :, (ii*K-K)+i,1)), squeeze(IRF.IRF_cs_ql(:,(ii*K-K)+i,1))],'k-')
                plot(horizon, -Plugin.IRFcum( i, : ),'k--', ...
                    horizon,-[InferenceMSW.MSWlboundcum(i, :)', InferenceMSW.MSWuboundcum(i,:)'],'k--')
                hold off
            else
                hold on
                plot(horizon, IRF.irfs_cs(:,(ii*K-K)+i ),'k-', ...
                    horizon,[squeeze(IRF.IRF_cs_qu( :, (ii*K-K)+i,1)), squeeze(IRF.IRF_cs_ql(:,(ii*K-K)+i,1))],'k-')
                plot(horizon, Plugin.IRFcum( i, : ),'k--', ...
                    horizon,[InferenceMSW.MSWlboundcum(i, :)', InferenceMSW.MSWuboundcum(i,:)'],'k--')
                hold off
            end
        else
            if negative(a)==1
                hold on
                plot(horizon, -IRF.irfs(:,(ii*K-K)+i ),'k-', ...
                    horizon,-[squeeze(IRF.IRF_qu( :, (ii*K-K)+i,1)), squeeze(IRF.IRF_ql(:,(ii*K-K)+i,1))],'k-')
                plot(horizon, -Plugin.IRF( i, : ),'k--', ...
                    horizon,-[InferenceMSW.MSWlbound(i, :)', InferenceMSW.MSWubound(i,:)'],'k--')
                hold off
            else
                hold on
                plot(horizon, IRF.irfs(:,(ii*K-K)+i ),'k-', ...
                    horizon,[squeeze(IRF.IRF_qu( :, (ii*K-K)+i,1)), squeeze(IRF.IRF_ql(:,(ii*K-K)+i,1))],'k-')
                plot(horizon, Plugin.IRF( i, : ),'k--', ...
                    horizon,[InferenceMSW.MSWlbound(i, :)', InferenceMSW.MSWubound(i,:)'],'k--')
                hold off
            end
            xlim([horizon(1),horizon(end)])
            title(strcat('$', shocks{ii},' \rightarrow ',variables{i},'$'),'Interpreter','latex')
        end
        a2 = a2 +1;
        xlim([horizon(1),horizon(end)])
        title(strcat('$', shocks{ii},' \rightarrow ',variables{i},'$'),'Interpreter','latex')
        if i==1
            ylabel(strcat('vs. SVAR-IV (', ivs_label(idx_IVs(a)),')'),'Interpreter','latex') 
        end
    end 
end
set(hFig1,'PaperPositionMode','auto')
set(hFig1, 'Position', [30 50  900 700])
hFig1 = tightfig(hFig1);  
print(hFig1,'output/figures/IRFsVsProxy', '-painters' ,'-dpdf') 







