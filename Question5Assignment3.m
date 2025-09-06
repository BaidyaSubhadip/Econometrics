clear all;
clc;
format compact; 
try
    [data, headers] = xlsread("time.xlsx");
catch ME
    error('Could not read "time.xlsx". Make sure the file is in the correct directory.');
end
% Assign variables from the data file
year = data(:,1);
y1 = data(:,2);
y2 = data(:,3);
y3 = data(:,4);
x2 = data(:,5);
x3 = data(:,6);
% --- Define common variables for all parts ---
Dt = double(year >= 1939 & year <= 1945);
Ct = 1 - Dt;
n = length(y1); % Number of observations
%% ----------------- Part 5a -----------------
fprintf('Question 5(a): Regress y1 on intercept, x2, x3\n');
fprintf('---------------------------------------------------\n');
X_a = [x2 x3];
stats_a = regstats(y1, X_a, 'linear', {'beta', 'tstat', 'rsquare', 'mse'});
beta_a = stats_a.beta;
se_a = stats_a.tstat.se;
tstat_a = stats_a.tstat.t;
rsquare_a = stats_a.rsquare;
stderr_reg_a = sqrt(stats_a.mse);
% Display the results
rowlabels_a = char('Intercept', 'x2', 'x3');
fprintf('%-10s %12s %12s %12s\n', 'Variable', 'Coefficient', 'Std. Error', 't-Statistic');
fprintf('---------------------------------------------------\n');
for i = 1:length(beta_a)
    fprintf('%-10s %12.4f %12.4f %12.4f\n', rowlabels_a(i,:), beta_a(i), se_a(i), tstat_a(i));
end
fprintf('---------------------------------------------------\n');
fprintf('R-squared: %.4f\n', rsquare_a);
fprintf('Standard Error of Regression: %.4f\n\n\n', stderr_reg_a);
%% ----------------- Part 5b -----------------
fprintf('Question 5(b): Regress y3 on intercept, x2, x3\n');
fprintf('---------------------------------------------------\n');
X_b = [x2 x3];
stats_b = regstats(y3, X_b, 'linear', {'beta', 'tstat', 'rsquare', 'mse'});
beta_b = stats_b.beta;
se_b = stats_b.tstat.se;
tstat_b = stats_b.tstat.t;
rsquare_b = stats_b.rsquare;
stderr_reg_b = sqrt(stats_b.mse);
% Display the results
rowlabels_b = char('Intercept', 'x2', 'x3');
fprintf('%-10s %12s %12s %12s\n', 'Variable', 'Coefficient', 'Std. Error', 't-Statistic');
fprintf('---------------------------------------------------\n');
for i = 1:length(beta_b)
    fprintf('%-10s %12.4f %12.4f %12.4f\n', rowlabels_b(i,:), beta_b(i), se_b(i), tstat_b(i));
end
fprintf('---------------------------------------------------\n');
fprintf('R-squared: %.4f\n', rsquare_b);
fprintf('Standard Error of Regression: %.4f\n\n\n', stderr_reg_b);
%% ----------------- Part 5c -----------------
fprintf('Question 5(c): Regress y1 on Ct, Dt, x2, x3 (No Intercept Model)\n');
fprintf('--------------------------------------------------------------------\n');
% Define Ct and Dt
Dt = double(year >= 1939 & year <= 1945);
Ct = 1 - Dt;

% Regressors for the no-intercept model
X_c = [Ct Dt x2 x3];

stats_c = regstats(y1, X_c, eye(4), {'beta', 'covb', 'tstat', 'rsquare', 'mse'});

% Extract values
beta_c = stats_c.beta;
se_c = stats_c.tstat.se;
tstat_c = stats_c.tstat.t;
rsquare_c = stats_c.rsquare;
stderr_reg_c = sqrt(stats_c.mse);
cov_matrix_c = stats_c.covb;

% Display the results
rowlabels_c = char('Ct', 'Dt', 'x2', 'x3');
fprintf('%-10s %12s %12s %12s\n', 'Variable', 'Coefficient', 'Std. Error', 't-Statistic');
fprintf('--------------------------------------------------------------------\n');
for i = 1:length(beta_c)
    fprintf('%-10s %12.4f %12.4f %12.4f\n', rowlabels_c(i,:), beta_c(i), se_c(i), tstat_c(i));
end
fprintf('--------------------------------------------------------------------\n');
fprintf('R-squared: %.4f\n', rsquare_c);
fprintf('Standard Error of Regression: %.4f\n\n', stderr_reg_c);
fprintf('Variance-Covariance Matrix:\n');
disp(cov_matrix_c);
fprintf('\nInterpretation Note for 5(c):\n');
fprintf(' - The coefficient on Ct (beta_11) is the estimated intercept for non-war years.\n');
fprintf(' - The coefficient on Dt (beta_21) is the estimated intercept for war years.\n');
fprintf(' - Including a separate intercept term would cause perfect multicollinearity because Ct + Dt = 1.\n\n\n');


%% ----------------- Part 5d -----------------
fprintf('Question 5(d): Calculate Var(delta_hat) where delta = beta_21 - beta_11\n');
fprintf('-----------------------------------------------------------------------\n');

% Var(δ̂) = Var(β̂_Dt - β̂_Ct) = Var(β̂_Dt) + Var(β̂_Ct) - 2*Cov(β̂_Dt, β̂_Ct)
% Note: In our matrix, Ct is the 1st variable and Dt is the 2nd.
var_ct = cov_matrix_c(1,1);
var_dt = cov_matrix_c(2,2);
cov_dt_ct = cov_matrix_c(1,2); % Or cov_matrix_c(2,1)

var_delta = var_dt + var_ct - 2 * cov_dt_ct;

fprintf('Var(β̂_Dt - β̂_Ct) = Var(β̂_Dt) + Var(β̂_Ct) - 2*Cov(β̂_Dt, β̂_Ct)\n');
fprintf('Var(δ̂) = %.4f + %.4f - 2*(%.4f) = %.4f\n\n\n', var_dt, var_ct, cov_dt_ct, var_delta);


%% ----------------- Part 5e -----------------
fprintf('Question 5(e): Regress y1 on intercept, Dt, x2, x3\n');
fprintf('--------------------------------------------------------------------\n');

% Regressors for the intercept model with a dummy
X_e = [Dt x2 x3];

stats_e = regstats(y1, X_e, 'linear', {'beta', 'tstat', 'rsquare', 'mse'});

% Extract values
beta_e = stats_e.beta;
se_e = stats_e.tstat.se;
tstat_e = stats_e.tstat.t;
rsquare_e = stats_e.rsquare;
stderr_reg_e = sqrt(stats_e.mse);

% Display the results
rowlabels_e = char('Intercept', 'Dt', 'x2', 'x3');
fprintf('%-10s %12s %12s %12s\n', 'Variable', 'Coefficient', 'Std. Error', 't-Statistic');
fprintf('--------------------------------------------------------------------\n');
for i = 1:length(beta_e)
    fprintf('%-10s %12.4f %12.4f %12.4f\n', rowlabels_e(i,:), beta_e(i), se_e(i), tstat_e(i));
end
fprintf('--------------------------------------------------------------------\n');
fprintf('R-squared: %.4f\n', rsquare_e);
fprintf('Standard Error of Regression: %.4f\n\n', stderr_reg_e);
fprintf('Interpretation Note for 5(e):\n');
fprintf(' - The `Intercept` is the baseline intercept for non-war years (when Dt=0).\n');
fprintf(' - The coefficient on `Dt` is the *difference* in the intercept between war years and non-war years.\n');
fprintf(' - The intercept for war years is (Intercept + Coef_Dt) = %.4f + %.4f = %.4f\n', beta_e(1), beta_e(2), beta_e(1)+beta_e(2));
fprintf(' - Note that this sum is identical to the coefficient on Dt in part (c), as it should be.\n');

%% ----------------- Part 5f -----------------
fprintf('Question 5(f): Regress y2 on intercept, Dt, x2, x3, xt4\n');
fprintf('--------------------------------------------------------------------\n');

% Create interaction variable xt4 = x3 * Dt
xt4 = x3 .* Dt;

% Regressors for the model
X_f = [Dt x2 x3 xt4];

stats_f = regstats(y2, X_f, 'linear', {'beta', 'tstat', 'rsquare', 'mse'});

% Extract values
beta_f = stats_f.beta;
se_f = stats_f.tstat.se;
rsquare_f = stats_f.rsquare;
stderr_reg_f = sqrt(stats_f.mse);

% Display results
rowlabels_f = char('Intercept', 'Dt', 'x2', 'x3', 'xt4');
fprintf('%-10s %12s %12s\n', 'Variable', 'Coefficient', 'Std. Error');
fprintf('--------------------------------------------------------------------\n');
for i = 1:length(beta_f)
    fprintf('%-10s %12.4f %12.4f\n', rowlabels_f(i,:), beta_f(i), se_f(i));
end
fprintf('--------------------------------------------------------------------\n');
fprintf('R-squared: %.4f\n', rsquare_f);
fprintf('Standard Error of Regression: %.4f\n\n\n', stderr_reg_f);


%% ----------------- Part 5g -----------------
fprintf('Question 5(g): Regress y2 on Ct, Dt, x2, xt4, xt5 & F-test\n');
fprintf('--------------------------------------------------------------------\n');

% Create interaction variables xt4 and xt5
xt4 = x3 .* Dt;
xt5 = x3 .* Ct;

% Regressors for the UNRESTRICTED model
X_g_unrestricted = [Ct Dt x2 xt4 xt5];
k_unrestricted = size(X_g_unrestricted, 2); % Number of parameters = 5

% Run the UNRESTRICTED regression (NO intercept since Ct and Dt are included)
stats_g = regstats(y2, X_g_unrestricted, eye(k_unrestricted), {'beta', 'tstat', 'rsquare', 'mse', 'covb', 'r'});

% Extract values
beta_g = stats_g.beta;
se_g = stats_g.tstat.se;
rsquare_g = stats_g.rsquare;
stderr_reg_g = sqrt(stats_g.mse);
cov_matrix_g = stats_g.covb;
SSE_U = sum(stats_g.r.^2); % Sum of Squared Errors for Unrestricted model

% Display results
rowlabels_g = char('Ct', 'Dt', 'x2', 'xt4', 'xt5');
fprintf('%-10s %12s %12s\n', 'Variable', 'Coefficient', 'Std. Error');
fprintf('--------------------------------------------------------------------\n');
for i = 1:length(beta_g)
    fprintf('%-10s %12.4f %12.4f\n', rowlabels_g(i,:), beta_g(i), se_g(i));
end
fprintf('--------------------------------------------------------------------\n');
fprintf('R-squared: %.4f\n', rsquare_g);
fprintf('Standard Error of Regression: %.4f\n\n', stderr_reg_g);
fprintf('Variance-Covariance Matrix:\n');
disp(cov_matrix_g);

fprintf('Explanation: Why x3 is not included in the model:\n');
fprintf('x3 is perfectly collinear with xt4 and xt5. Specifically, xt4 + xt5 = (x3.*Dt) + (x3.*Ct) = x3.*(Dt+Ct) = x3.*1 = x3.\n');
fprintf('Including all three would violate the full-rank assumption of OLS (dummy variable trap for interactions).\n\n');

% --- F-Test for H0: beta_4 = beta_5 = 0 ---
fprintf('--- F-Test for H0: beta_4 = 0 and beta_5 = 0 ---\n');

% Run the RESTRICTED model (y2 on Ct, Dt, x2)
X_g_restricted = [Ct Dt x2];
k_restricted = size(X_g_restricted, 2); % Number of parameters = 3
stats_g_restricted = regstats(y2, X_g_restricted, eye(k_restricted), {'r'});
SSE_R = sum(stats_g_restricted.r.^2); % Sum of Squared Errors for Restricted model

q = k_unrestricted - k_restricted; % Number of restrictions = 2

% Calculate F-statistic
F_statistic = ((SSE_R - SSE_U) / q) / (SSE_U / (n - k_unrestricted));
p_value = 1 - fcdf(F_statistic, q, n - k_unrestricted);
alpha = 0.05;
F_critical = finv(1 - alpha, q, n - k_unrestricted);

fprintf('Null Hypothesis (H0): The coefficients for xt4 and xt5 are both zero.\n');
fprintf('SSE (Restricted) = %.4f\n', SSE_R);
fprintf('SSE (Unrestricted) = %.4f\n', SSE_U);
fprintf('F-statistic = %.4f\n', F_statistic);
fprintf('p-value = %.4f\n', p_value);
fprintf('Critical F-value at alpha=%.2f (df=%d, %d) = %.4f\n', alpha, q, n-k_unrestricted, F_critical);

if F_statistic > F_critical
    fprintf('Conclusion: Since F-statistic > F-critical (%.4f > %.4f), we reject the null hypothesis.\n', F_statistic, F_critical);
else
    fprintf('Conclusion: Since F-statistic <= F-critical (%.4f <= %.4f), we fail to reject the null hypothesis.\n', F_statistic, F_critical);
end
fprintf('There is statistically significant evidence that at least one of the interaction terms (xt4, xt5) is non-zero.\n\n\n');


%% ----------------- Part 5h -----------------
fprintf('Question 5(h): Regress y3 on intercept, Dt, x2, x3, xt4, xt6\n');
fprintf('--------------------------------------------------------------------\n');

% Create the new interaction variable xt6 = x2 * Dt
xt6 = x2 .* Dt;
% We already have xt4 = x3 .* Dt from part (f)

% Regressors for the model
X_h = [Dt x2 x3 xt4 xt6];

% Use 'linear' to automatically include an intercept
stats_h = regstats(y3, X_h, 'linear', {'beta', 'tstat', 'rsquare', 'mse'});

% Extract values
beta_h = stats_h.beta;
se_h = stats_h.tstat.se;
rsquare_h = stats_h.rsquare;
stderr_reg_h = sqrt(stats_h.mse);

% Display results
rowlabels_h = char('Intercept', 'Dt', 'x2', 'x3', 'xt4', 'xt6');
fprintf('%-10s %12s %12s\n', 'Variable', 'Coefficient', 'Std. Error');
fprintf('--------------------------------------------------------------------\n');
for i = 1:length(beta_h)
    fprintf('%-10s %12.4f %12.4f\n', rowlabels_h(i,:), beta_h(i), se_h(i));
end
fprintf('--------------------------------------------------------------------\n');
fprintf('R-squared: %.4f\n', rsquare_h);
fprintf('Standard Error of Regression: %.4f\n\n', stderr_reg_h);