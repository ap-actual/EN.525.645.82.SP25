% EN.525.645.82.SP25 Homework 7
% Brian Caskey 

%% 7.1.1
% Solve a 3 x 3 system of simultaneous linear equations that you invent, 
% using Cramer’s rule by hand.  If you know how to invert a matrix by hand,
% or if you have the software that is capable of matrix inversion, use that
% technique as well. Then, check your results, by hand.

clear,clc

% Coefficient matrices
A = [2 1 -1;
    -3 -1 2;
    -2 1 2];

b = [8; -11; -3];

% determinant of A
detA = det(A);

% replace each column of A with b to get Ax, Ay, Az
Ax = A; Ax(:,1) = b;
Ay = A; Ay(:,2) = b;
Az = A; Az(:,3) = b;

% compute determinants
detAx = det(Ax);
detAy = det(Ay);
detAz = det(Az);

% apply Cramer's Rule
x = detAx / detA;
y = detAy / detA;
z = detAz / detA;

% display 
fprintf('x = %f\n', x);
fprintf('y = %f\n', y);
fprintf('z = %f\n', z);


%% 7.2.1 
% Invent a system of three equations in two unknowns. Graph the equations 
% and find the least squares solution. Plot the solution on your graph.

% my equations are: 
%     x + 2y = 5
%    2x -  y = 3
%    -x +  y = 1

% Define the system
A = [ 1  2; 
      2 -1; 
     -1  1];

b = [5; 
     3; 
     1];

% solve using least squares
x_ls = (A' * A)^-1 * A' * b;  

% vals to plot
x = linspace(-2, 5, 100);

% solve for y vals for plot (after re-arranging equations above)
y1 = (5 - x)/2;
y2 = (2*x - 3);
y3 = (1 + x);

% plot
figure;
plot(x, y1, 'LineWidth', 2); hold on;
plot(x, y2, 'LineWidth', 2);
plot(x, y3, 'LineWidth', 2);

% plot least squares solution
plot(x_ls(1), x_ls(2), 'kx', 'MarkerSize', 10);

legend('x + 2y = 5', '2x - y = 3', '-x + y = 1', 'Least Squares Solution');
xlabel('x'); ylabel('y');
grid on;
axis equal;

fprintf('X = %f\n', x_ls(1));
fprintf('Y = %f\n', x_ls(2));

%% 7.2.2 
% Using the example in the presentation, compute the coordinates x1, y1 by
% differentiating the expression for the sum of the squares of the 
% distances to each of the three lines.

A = [ 3 -1; 
     -1 3];
b = [4;
     2];

x_lsq = A^-1 * b

%% 7.3.1 
% Using random numbers, and setting the deterministic error to zero, 
% simulate equations 3.1 through 3.6 of the text.

clear,clc

% reset random number generator for repeatable results
rng('default');

% problem assumptions
i      = 1:100;             % steps to run
N      = numel(i);
eta_di = zeros(1,N); % deterministic error to 0

% truth data
x_a = 3.0 * ones(1, N);

% compute random measurement variation
mean_eta_i_1 = 0;
std_eta_i_1  = 1;
eta_i_1  = mean_eta_i_1 + std_eta_i_1 .* randn(1,N);

% eq 3.1
x_i = x_a + eta_di + eta_i_1;

% eq 3.3
% The best estimate of position in this case is the average of the measurements
x_hat_N = 1/N * sum(x_i);

% where the standard deviation of the error in the estimated position 
delta_x_hat_N = x_hat_N - x_a; % eq 3.4

% is the standard deviation of the random error sigma_eta divided by the
% square root of the number of measurements
sigma_delta_x_hat_N = std_eta_i_1 / (sqrt(N));

% hence, the sequence of measurements can be processed to obtain and
% ever-improving estiamte of the vehicles position. To avoid storing all
% past measurements, the Nth estimate can be written in the recursive form:

x_hat_iN = zeros(1,N);
for iN = 1:N-1
    x_hat_iN(iN+1) = ((iN - 1)/iN) * x_hat_iN(iN) + x_i(iN)/iN;
end


% plot!
figure(); hold on; grid on;
plot(x_a, 'b.', 'MarkerSize', 10); lgnd(1) = "Truth Measurements";
plot(x_i, 'r.', 'MarkerSize', 7); lgnd(end+1) = "Position Measurements";

c="#A2142F";
lw = 1.5;
yline(x_hat_N, 'Color', c, 'LineWidth', lw); lgnd(end+1) = "Measurement Average";
yline(x_hat_N+sigma_delta_x_hat_N, '--', 'Color', c, 'LineWidth', lw); lgnd(end+1) = "Measurement Average 1-σ";
yline(x_hat_N+-1*sigma_delta_x_hat_N, '--', 'Color', c, 'LineWidth', lw); lgnd(end+1) = "";

plot(x_hat_iN, ':', 'Color', c, 'LineWidth', lw); lgnd(end+1) = "Recursive form";

legend(lgnd)

%% 7.3.2
% Extend this to simulation of equations 3.7 to 3.11. Make the simulated 
% errors in one set of measurements significantly larger than those in the 
% other, thus illustrating that inclusion of measurements with low weights 
% is in principle helpful, but in practice is not always worth doing

% compute random measurement variation
mean_eta_i_2 = 0;
std_eta_i_2  = 20;
eta_i_2  = mean_eta_i_2 + std_eta_i_2 .* randn(1,N);

% eq 3.1
x_i_2 = x_a + eta_di + eta_i_2;

% eq 3.3
% The best estimate of position in this case is the average of the measurements
x_hat_N2 = 1/N * sum(x_i_2);

% compute weights - eq 3.8
D_N = std_eta_i_1^2 + std_eta_i_2^2;
lambda_1 = std_eta_i_2^2 / D_N;
lambda_2 = std_eta_i_1^2 / D_N;

% compute best estimate
x_bar_N   = lambda_1 * x_hat_N + lambda_2*x_hat_N2;
x_bar_var =  std_eta_i_1^2* std_eta_i_2^2 / D_N;

% plot!
figure(); hold on; grid on; clear lgnd
plot(x_a, 'b.', 'MarkerSize', 10); lgnd(1) = "Truth";
plot(x_i, 'r.', 'MarkerSize', 7); lgnd(end+1) = "Position Measurements set 1";
plot(x_i_2, 'rx', 'MarkerSize', 7); lgnd(end+1) = "Position Measurements set 2";

c="#A2142F";
lw = 1.5;
yline(x_bar_N, 'Color', c, 'LineWidth', lw); lgnd(end+1) = "Best Estimate via WLSq";
yline(x_hat_N+sqrt(x_bar_var), '--', 'Color', c, 'LineWidth', lw); lgnd(end+1) = "1-σ";
yline(x_hat_N+-1*sqrt(x_bar_var), '--', 'Color', c, 'LineWidth', lw); lgnd(end+1) = "";

legend(lgnd)
 

