%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prevalence estimation with non-random testing
% Benatia, Godefroy and Lewis (2020)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DATA is the complete dataset
% nt is nb of daily tests/population
% nt_2 is nb of daily tests
% lag_nt is lagged nb of daily tests/population
% st is nb of daily positive tests
% yt is the dependent variable (diff of daily log positive case rates)

% Select observation such that:
b = find(DATA.nt>0 & DATA.st>0 & DATA.lag_nt>0 & DATA.lag_st>0);

% Keep only obs satisfying above criterion:
DATA_r.nt = DATA.nt(b);
DATA_r.nt_2 = DATA.nt_2(b);
DATA_r.st = DATA.st(b);
DATA_r.lag_nt = DATA.lag_nt(b);
DATA_r.lag_st = DATA.lag_st(b);
DATA_r.yt = DATA.yt(b);
DATA_r.date = DATA.date(b);
DATA_r.pop = DATA.pop(b);
DATA_r.state = DATA.state(b);

% Estimation procedure (needs criterion function log_lik2.m, see below)

options = optimoptions('fminunc','SpecifyObjectiveGradient',false,... % Options
    'Algorithm','quasi-newton',...
    'Display','none',...
    'OptimalityTolerance',1e-10,...
    'StepTolerance',1e-10,...
    'MaxFunctionEvaluations',10000,...
    'CheckGradients',false);
PARAM_init = [3,3,3,-1000,1]; %Set initial parameter values for optimization routine

[beta,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] =  fminunc(@(b0) loglik_2(b0,DATA_r),PARAM_init,options);
COV_beta = inv(HESSIAN);

disp( ' %%%%%%%%%%%%%%%%%%%% Results')
disp([' alpha1  ',num2str(beta(1)),' (',num2str(sqrt(COV_beta(1,1))),')']);
disp([' alpha2  ',num2str(beta(2)),' (',num2str(sqrt(COV_beta(2,2))),')']);
disp([' alpha3  ',num2str(beta(3)),' (',num2str(sqrt(COV_beta(3,3))),')']);
disp([' beta  ',num2str(beta(4)),' (',num2str(sqrt(COV_beta(4,4))),')']);
disp([' sigma_u  ',num2str(beta(5)),' (',num2str(sqrt(COV_beta(5,5))),')']);
disp(['observations ',num2str(length(DATA_r.nt))])
disp( ' %%%%%%%%%%%%%%%%%%%% ')

beta_hat = beta;

% Estimate prevalence from above estimates
DATA_r.prevalence = exp( log((DATA_r.st./DATA_r.nt_2))  ...
    +beta_hat(1)*( exp(beta_hat(4)) - exp(beta_hat(4)*DATA_r.nt) ) ...
    +beta_hat(2)*( exp(2*beta_hat(4)) - exp(2*beta_hat(4)*DATA_r.nt) ) ...
    +beta_hat(3)*( exp(3*beta_hat(4)) - exp(3*beta_hat(4)*DATA_r.nt) ) ...
    );

% Compute confidence interval
ds = ones(length(DATA_r.pop),1);
dn = - (beta(1)*1*beta(4)*exp(1*beta(4)*DATA_r.nt)./DATA_r.pop)...
    - (beta(2)*1*beta(4)*exp(2*beta(4)*DATA_r.nt)./DATA_r.pop)...
    - (beta(3)*1*beta(4)*exp(3*beta(4)*DATA_r.nt)./DATA_r.pop);
dalpha1 = exp(1*beta(4)) - exp(1*beta(4)*DATA_r.nt);
dalpha2 = exp(2*beta(4)) - exp(2*beta(4)*DATA_r.nt);
dalpha3 = exp(3*beta(4)) - exp(3*beta(4)*DATA_r.nt);
dbeta = beta(1)*1*( exp(1*beta(4)) - DATA_r.nt.*exp(1*beta(4)*DATA_r.nt)   )...
    +beta(2)*2*( exp(2*beta(4)) - DATA_r.nt.*exp(2*beta(4)*DATA_r.nt)   )...
    +beta(3)*3*( exp(3*beta(4)) - DATA_r.nt.*exp(3*beta(4)*DATA_r.nt)   );

H = [cov([DATA_r.nt,log(DATA_r.st)]),zeros(2,4);...
    zeros(4,2),COV_beta(1:4,1:4)];

dlogS = [ds,dn,dalpha1,dalpha2,dalpha3,dbeta];
H = sqrt(diag(dlogS*H*dlogS'));

logS =  log((DATA_r.st./DATA_r.nt_2))  ...
    +beta_hat(1)*( exp(beta_hat(4)) - exp(beta_hat(4)*DATA_r.nt) ) ...
    +beta_hat(2)*( exp(2*beta_hat(4)) - exp(2*beta_hat(4)*DATA_r.nt) ) ...
    +beta_hat(3)*( exp(3*beta_hat(4)) - exp(3*beta_hat(4)*DATA_r.nt) ) ...
    ;

logS_hi = logS+1.96*H;
logS_lo = logS-1.96*H;

DATA_r.prevalence_hi = exp(logS_hi);
DATA_r.prevalence_lo = exp(logS_lo);
