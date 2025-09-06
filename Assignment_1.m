clear;
clc;
[data, headers] = xlsread('ProdFuncData.xlsx');
valueadded = data(:,2);
labor = data(:,3);
capital = data(:,4);
% ---Part 5(a)---%
y = log(valueadded);           % Take log of output (log Y)
loglabor = log(labor);         % log(L)
logcapital = log(capital);     % log(K)
x = [loglabor logcapital];
whichstats = {'beta','covb','yhat','r','rsquare','mse','tstat','fstat'};
model = regstats(y, x, 'linear', whichstats);
betaM = model.beta;
seM = model.tstat.se;
tstatM = model.tstat.t;
covM = model.covb;
rsquare = model.rsquare;
rowlabels = char('Intercept','LogLabor','LogCapital');
fprintf('              Model for 5.a               \n')
fprintf(' ___________________________________\n')
fprintf(' %-12s  %8s  %7s  %6s\n','Variable','Est','se','tstat')
fprintf(' %-12s  %8s  %7s  %6s\n','__________','________','______','______')
for i = 1:length(betaM)
    fprintf(' %-12s  %8.4f  %7.4f  %6.2f\n', ...
        rowlabels(i,:), betaM(i), seM(i), tstatM(i))
end
fprintf('\nCovariance Matrix of beta (Cov[b]):\n')
disp(covM)
fprintf('R-squared: %.4f\n', rsquare)
% ---Part 5.b----%
R = [0 1 1];  % Tests beta_L + beta_K = 1
q = 1;
J = 1;
Rbeta_minus_q = R * betaM - q;
RCovR = R * covM * R';
Fstat = (1/J) * (Rbeta_minus_q') * (inv(RCovR)) * (Rbeta_minus_q);
fprintf('The F-statistics for the hypothesis for the question 5b:')
disp(Fstat)
% ---Part 5(c)---%
logL2 = 0.5 * (loglabor.^2);
logK2 = 0.5 * (logcapital.^2);
logLK = loglabor .* logcapital;
X_translog = [loglabor, logcapital, logL2, logK2, logLK];
model_translog = regstats(y, X_translog, 'linear', whichstats);
betaT = model_translog.beta;
seT = model_translog.tstat.se;
tstatT = model_translog.tstat.t;
covT = model_translog.covb;
rsquareT = model_translog.rsquare;
rowlabelsT = char('Intercept', 'LogLabor', 'LogCapital', ...
                  '0.5*LogLabor^2', '0.5*LogCapital^2', 'LogLabor*LogCapital');
fprintf('\n              Translog Model for 5.c               \n')
fprintf(' _________________________________________________\n')
fprintf(' %-20s  %8s  %8s  %7s\n','Variable','Est','se','tstat')
fprintf(' %-20s  %8s  %8s  %7s\n','____________________','________','________','_______')

for i = 1:length(betaT)
    fprintf(' %-20s  %8.4f  %8.4f  %7.2f\n', ...
        rowlabelsT(i,:), betaT(i), seT(i), tstatT(i))
end
fprintf('\nCovariance Matrix of beta (Translog model):\n')
disp(covT)
fprintf('R-squared (Translog): %.4f\n', rsquareT)
% ---Part 5(d)---%
R = [0 0 0 1 0 0; 
     0 0 0 0 1 0; 
     0 0 0 0 0 1]; 
q = [0; 0; 0];        % H0: beta4 = beta5 = beta6 = 0
J = 3;        % Number of restrictions = 3
Rbeta_minus_q = R * betaT - q;
RCovR = R * covT * R';
Fstat_d = (1/J) * (Rbeta_minus_q') * (inv(RCovR)) * (Rbeta_minus_q);
fprintf('\n--- Part 5.d ---\n')
fprintf('The F-statistics for the hypothesis for the question 5d: %.4f\n', Fstat_d)
% ---Part 5(e)---%
R = [0 1 1 0 0 0;
     0 0 0 1 1 2];
q = [1; 0];
J = 2;
Rbeta_minus_q = R * betaT - q;
RCovR = R * covT * R';
Fstat_e = (1/J) * (Rbeta_minus_q') * (inv(RCovR)) * (Rbeta_minus_q);
fprintf('\n--- Part 5.e ---\n')
fprintf('The F-statistics for the hypothesis for the question 5e: %.4f\n', Fstat_e)


