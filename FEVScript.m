
%--------------------------------------------------------------------------
% This example from Bayesian Ideas and Data Analysis by Wesley Johnson
% etal.
%--------------------------------------------------------------------------

clear all
clc
[data, headers] = xlsread('/Users/arshadrahman/Library/CloudStorage/Dropbox/1-WorkFiles/1 IIT-Kanpur/Courses/ECO251 Summer 2025/Practice1/FEV Wesley.xlsx',1,'A2:C343');
% The study is for individuals aged 10 to 18. Hence the last 3 observations
% have been left out.

age   = data(:,1);
fev   = data(:,2);
smoke = data(:,3);

y        = fev;    % forced expiratory volume
agesqr   = age.^2;
ageIntSmoke = age.*smoke;

%% Data Summary
disp('__________ Smoke __________')
tabulate(smoke)

ageSmoke     = age(find(smoke==1));
ageNonSmoke  = age(find(smoke~=1));

fevSmoke     = fev(find(smoke==1));
fevNonSmoke  = fev(find(smoke~=1));

meanage = [mean(ageSmoke); mean(ageNonSmoke)];
meanfev = [mean(fevSmoke); mean(fevNonSmoke)];

stdage = [std(ageSmoke); std(ageNonSmoke)];
stdfev = [std(fevSmoke); std(fevNonSmoke)];

%----- Part 1(a) -------

rowlabels = char('Smoke = 1','Smoke = 0');
% Data Summary

        fprintf('                  FEV                 Age \n        ')
        fprintf('    ________________   _________________  \n')
        fprintf('             mean      std      mean       std     \n')
        for i = 1:2
         fprintf(' %-8s  % 4.2f  % 8.2f  % 8.2f  % 8.2f  \n', rowlabels(i,:), meanfev(i), ...
             stdfev(i), meanage(i), stdage(i))
        end

%% Regression

x1 = [age smoke];
x2 = [age agesqr smoke];
x3 = [age agesqr smoke ageIntSmoke];

whichstats ={'beta','covb','yhat','r','rsquare','mse','tstat','fstat'};

stats1 = regstats(y,x1,'linear',whichstats);
stats2 = regstats(y,x2,'linear', whichstats);
stats3 = regstats(y,x3,'linear',whichstats);

betaM1 = stats1.beta;
betaM2 = stats2.beta;
betaM3 = stats3.beta;

seM1   = stats1.tstat.se;
seM2   = stats2.tstat.se;
seM3   = stats3.tstat.se;

tstatM1   = stats1.tstat.t;
tstatM2   = stats2.tstat.t;
tstatM3   = stats3.tstat.t;

%----- Part 1(b) -------

rowlabels1 = char('Intercept','age','smoke');
     fprintf('              Model M1              \n')
     fprintf(' ___________________________________\n')
     fprintf('             Est     se    tstat    \n')
     fprintf('            ______  _____  ______   \n')

     for i = 1:length(betaM1)
         fprintf(' %-8s  % 4.2f  % 4.2f  % 6.2f  \n', ...
           rowlabels1(i,:), betaM1(i), seM1(i), tstatM1(i))
     end

%% ----- Part 1(c) -------

rowlabels2 = char('Intercept','age','agesqr','smoke');
     fprintf('              Model M2              \n')
     fprintf(' ___________________________________\n')
     fprintf('             Est     se    tstat    \n')
     fprintf('            ______  _____  ______   \n')

     for i = 1:length(betaM2)
         fprintf(' %-8s  % 4.2f  % 4.2f  % 6.2f  \n', ...
           rowlabels2(i,:), betaM2(i), seM2(i), tstatM2(i))
     end
     
%% ----- Part 1(d) -------

rowlabels3 = char('Intercept','age','agesqr','smoke','age*smoke');
     fprintf('              Model M3              \n')
     fprintf(' ___________________________________\n')
     fprintf('             Est     se    tstat    \n')
     fprintf('            ______  _____  ______   \n')

     for i = 1:length(betaM3)
         fprintf(' %-8s  % 4.2f  % 4.2f  % 6.2f  \n', ...
           rowlabels3(i,:), betaM3(i), seM3(i), tstatM3(i))
     end
     
rsq1 = stats1.rsquare   
rsq2 = stats2.rsquare     
rsq3 = stats3.rsquare     



%% ----- Part 1(e) -------

x2plot = min(age):0.01:max(age);
x3plot = x2plot.^2;
m      = size(x2plot,2);

x1plot        = ones(m,1);
x4smokeplot   = ones(m,1);
x5smokeplot   = x2plot'.*x4smokeplot;

x4nonsmokeplot   = zeros(m,1);
x5nonsmokeplot = x2plot'.*x4nonsmokeplot;

xsmokeplot    = [x1plot x2plot' x3plot' x4smokeplot x5smokeplot];
xnonsmokeplot = [x1plot x2plot' x3plot' x4nonsmokeplot x5nonsmokeplot];

ysmokeplot    = zeros(m,1);
ynonsmokeplot = zeros(m,1);

for i = 1:m
    ysmokeplot(i)    = xsmokeplot(i,:)*stats3.beta;
    ynonsmokeplot(i) = xnonsmokeplot(i,:)*stats3.beta;
end

plot(x2plot,ysmokeplot', 'color', [0.6,0.6,0.6], 'LineWidth',0.9,'LineStyle','-')
xlim([9,20])
text(max(x2plot)+0.25,max(ysmokeplot),'smoke=1')
hold on
plot(x2plot,ynonsmokeplot,'color', 'k', 'LineWidth',0.9,'LineStyle','-')
text(max(x2plot)+0.25,max(ynonsmokeplot),'smoke=0')
title('FEV data: Conditional Mean Functions')
xlabel('Age (year)','FontSize', 9)
ylabel('FEV','FontSize', 9)

hold on
scatter(age,y,[],[0.5,0.5,0.5])


%----- Part 1(f) -------

F = ((stats3.rsquare - stats1.rsquare)/2)/((1 - stats3.rsquare)/(size(fev,1)-5))   
% F = 6.7122, Ftable = finv(0.95,2,342 - 5) = 3.0225  


