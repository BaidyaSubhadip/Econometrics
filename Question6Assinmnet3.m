clear all;
clc;
format compact;
try
    data = xlsread("Advertising.xlsx");
catch ME
    error('Could not read "Advertising.xlsx". Make sure the file is in the correct directory.');
end
y = data(:,1); 
x = data(:,2); 
n = length(y);
%% ----------------- Part (a) -----------------
fprintf('Question 6(a): Fit simple linear regression model\n');
fprintf('---------------------------------------------------\n');
stats_a = regstats(y, x, 'linear', {'beta', 'tstat', 'r', 'mse'});

beta_a = stats_a.beta;
se_a = stats_a.tstat.se;
residuals = stats_a.r;
rowlabels_a = char('Intercept', 'x');
fprintf('%-10s %12s %12s\n', 'Variable', 'Coefficient', 'Std. Error');
fprintf('---------------------------------------------------\n');
for i = 1:length(beta_a)
    fprintf('%-10s %12.4f %12.4f\n', rowlabels_a(i,:), beta_a(i), se_a(i));
end
fprintf('---------------------------------------------------\n');
fprintf('Standard Error of Regression: %.4f\n\n', sqrt(stats_a.mse));
fprintf('Residuals (all %d many):\n', n);
disp(residuals(:));
fprintf('\n\n');

%% ----------------- Part (b) -----------------
fprintf('Question 6(b): Plot residuals against time\n');
fprintf('---------------------------------------------------\n');
time_index = (1:n)';
figure;
scatter(time_index, residuals, 'b', 'filled');
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
hold off;
title('Residuals vs. Time');
xlabel('Time (Month Index)');
ylabel('Residuals');
grid on;
%% ----------------- Part (c) -----------------
fprintf('Question 6(c): Formal test for positive correlation (Durbin-Watson)\n');
fprintf('--------------------------------------------------------------------\n');

% Manual calculation of the Durbin-Watson statistic
% d = sum( (e_t - e_{t-1})^2 ) / sum( e_t^2 )
numerator_dw = sum(diff(residuals).^2);
denominator_dw = sum(residuals.^2);
d_statistic = numerator_dw / denominator_dw;
fprintf('Manually Calculated Durbin-Watson Statistic (d) = %.4f\n\n', d_statistic);

