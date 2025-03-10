function [sigma] = ComputeSize(Phi_1)
% Extract data
x = 1:length(Phi_1{2}(65,:)); % Define x-axis values
y_1 = abs(Phi_1{1}(65,:)).^2;  % Compute squared magnitude
y_2 = abs(Phi_1{2}(65,:)).^2;  % Compute squared magnitude

% Define Gaussian model
gaussEqn = 'a*exp(-((x-b)^2)/(2*c^2))'; 

% Initial guesses: 
startPoints = [max(y_1), mean(x), std(x)]; 

% Fit the data
fitResult_1 = fit(x(:), y_1(:), gaussEqn, 'Start', startPoints);
fitResult_2 = fit(x(:), y_2(:), gaussEqn, 'Start', startPoints);

% Extract sigma (c in the equation)
sigma_1 = fitResult_1.c;
sigma_2 = fitResult_2.c;

sigma = max(sigma_1,sigma_2);
end