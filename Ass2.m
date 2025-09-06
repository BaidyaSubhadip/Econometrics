% Observed data
y = [5,0,1,1,0,3,2,3,4,1];
n = length(y);
sum_y = sum(y);


denominator = prod(1 ./ factorial(y));

% Theta grids
theta1 = 0.8:0.05:3.8;
theta2 = 0.5:0.05:3.5;

likelihood = denominator .* exp(-n.*theta1) .* theta1.^(sum_y);

log_likelihood = -n.*theta2 + sum_y.*log(theta2) - sum(log(factorial(y))) + 30;

% Plot
figure;
yyaxis left
plot(theta1, likelihood, 'b-', 'LineWidth', 2);
ylabel('Likelihood L(\theta|y)');

yyaxis right
plot(theta2, log_likelihood, 'r--', 'LineWidth', 2);
ylabel('log-likelihood + 30');

xlabel('\theta');
title('Likelihood and Log-Likelihood (Full Formula)');
legend('Likelihood','Log-likelihood + 30');
grid on;
