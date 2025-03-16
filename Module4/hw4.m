% EN.525.645.82.SP25 Homework 4
% Written by Brian Caskey on 03/15/2025

clear, clc

%% 4.1.1 
% Use Newton’s method to find the value of x for which the curve 
% y = 2x + 3x^2 is equal to 0.
%
% note that y' = 2 + 6x

% initial guess
x = nan(1,1000);
x(1) = 5;

% convergence criteria
e = 0.001;

% run loop (stop at 1000 loops)
for n = 1:1000

    % compute f(x_n) and f'(x_n)
    f_x_n  = 2*x(n) + 3*x(n)^2;
    fp_x_n = 2 + 6*x(n);

    % now compute x_n+1 via the equation from 
    % https://en.wikipedia.org/wiki/Newton%27s_method
    %
    % x_n+1 = x_n - f(x_n) / f'(x_n)
    x(n+1) = x(n) - f_x_n / fp_x_n;

    % display iteration
    fprintf("x_n = %f\nx(n+1) = %f\n\n", x(n), x(n+1));

    % if the convergence criteria is met, finish loop
    if abs(x(n) - x(n+1)) < e
        break
    end

end


% ---- Run plot check for fun -----
figure(); grid on; hold on;

% set up x axis data
xPlt = -10:0.1:10;

% plot actual function
fnEvl = 2.*xPlt + 3.*xPlt.^2;
plot(xPlt, fnEvl, 'k', 'LineWidth', 1.3);

yline(0)

% now plot lines
c = colormap(turbo(numel(find(~isnan(x)))));
for i = 1:numel(find(~isnan(x)))
    
    % plot point on function
    scatter(x(i), 2*x(i) + 3*x(i)^2, 'x', 'MarkerEdgeColor', c(i,:)); 

    % plot tangent line
    m = 2 + 6*x(i);
    plot(xPlt,  m*(xPlt - x(i)) + 2*x(i) + 3*x(i)^2, '--', 'Color', c(i,:));
    
    xline(x(i+1), ':', 'Color', c(i,:), 'LineWidth', 1.5)
end


%% 4.1.2
% For M = E - e sin E, solve for E when M = pi/2 and e = .3
clear, clc

% From the lecture notes: 
% To solve for E, solve for the roots of the equation `f(E) = E - e sinE - M`
% - Since sin(E) can be expressed as a polynomial by using Taylor series, we know that f(E) has a solution
% - More importantly, this solution can be found using iterative techniques
%
% 1. write as f(x) = 0
% 2. guess an initial value x = x0
% 3. compute the derivative f'(x0)
% 4. solve iteratively for the value of x_n that makes f(x_n) = 0

% constants from problem
M = pi/2;
e = 0.3;

% convergence criteria
eps = 0.001;

% initial guess
x = nan(1,1000);
x(1) = 5;

% run loop (stop at 1000 loops)
for n = 1:1000

    % compute f(x_n) and f'(x_n)
    f_x_n  = x(n) - e * sin(x(n)) - M;
    fp_x_n = 1 - e * cos(x(n));

    % now compute x_n+1 via the equation from 
    % https://en.wikipedia.org/wiki/Newton%27s_method
    %
    % x_n+1 = x_n - f(x_n) / f'(x_n)
    x(n+1) = x(n) - f_x_n / fp_x_n;

    % display iteration
    fprintf("x_n = %f\nx(n+1) = %f\n\n", x(n), x(n+1));

    % if the convergence criteria is met, finish loop
    if abs(x(n) - x(n+1)) < eps
        break
    end

end


% ---- Run plot check for fun -----
figure(); grid on; hold on;

% set up x axis data
xPlt = -10:0.1:10;

% plot actual function
fnEvl = xPlt - e * sin(xPlt) - M;
plot(xPlt, fnEvl, 'r', 'LineWidth', 1.3);

yline(0)

% now plot lines
c = colormap(winter(numel(find(~isnan(x)))));
for i = 1:numel(find(~isnan(x)))
    
    % plot point on function
    scatter(x(i), x(i) - e * sin(x(i)) - M, 'x', 'MarkerEdgeColor', c(i,:)); 

    % plot tangent line
    m = 1 - e * cos(x(i));
    plot(xPlt,  m*(xPlt - x(i)) + x(i) - e * sin(x(i)) - M, '--', 'Color', c(i,:));

    xline(x(i+1), ':', 'Color', c(i,:), 'LineWidth', 1.5)
end


%% 4.2.2 
% Define two circles in a two-dimensional plane that overlap. Specify the 
% centers of the two circles as points x1, y1 and x2, y2. Use Newton’s 
% method to find the two points x01, y01 and x02, y02 at which the circles 
% intersect. This is how GPS works. Specifically, the locations of the GPS 
% satellites are known. A GPS receiver measures the distance to each 
% satellite, thus defining a circle of position on the surface of the 
% earth that corresponds to each satellite. The points at which the circles 
% from multiple satellites intersect are computed, thus defining one’s 
% location on the surface of the earth.
clear, clc

% circle equations
%   Circle 1: F1 = (x - x1)² + (y - y1)² - r1² = 0
%   Circle 2: F2 = (x - x2)² + (y - y2)² - r2² = 0
%
% The jacobian of the systen is
%   J = [dF1/dx  dF1/dy
%        dF2/dx  dF2/dy]
%
% where (using the chain rule) we can find
%   dF1/dx = 2*(x-x1)
%   dF1/dy = 2*(y-y1)
%   dF2/dx = 2*(x-x2)
%   dF2/dy = 2*(y/y2)
%
% The equations from the wiki page below are then followed
% https://en.wikipedia.org/wiki/Newton%27s_method#Multidimensional_formulations 

% circle 1
x1 = 5;
y1 = 10;
r1 = 9;

% circle 2
x2 = -10; 
y2 = -8;
r2 = 15;

% initial guess
X = nan(2,1000);
X(:,1) = [x1; 0];

% convergence criteria
eps = 0.001;

% run loop (stop at 1000 loops)
for n = 1:1000

    % compute circle equations
    F1 = (X(1,n) - x1)^2 + (X(2,n) - y1)^2 - r1^2;
    F2 = (X(1,n) - x2)^2 + (X(2,n) - y2)^2 - r2^2;

    % Jacobian matrix
    J = [2*(X(1,n) - x1), 2*(X(2,n) - y1);
         2*(X(1,n) - x2), 2*(X(2,n) - y2)];

    % Function vector
    F = [F1; F2];

    % Newton & update step: solve J * delta = -F
    delta = -J \ F;
    X(1,n+1) = X(1,n) + delta(1);
    X(2,n+1) = X(2,n) + delta(2);

    % display iteration
    fprintf("X(n) = %f %f \nX(n+1) = %f %f\n\n", X(1,n), X(2,n), X(1,n+1), X(2,n+1));

    % if convergence criteria is met, break
    if norm(delta) < eps
        break;
    end
 
end

% -- Now plot to confirm
th = 0:pi/50:2*pi;

xCir1 = r1 * cos(th) + x1; yCir1 = r1 * sin(th) + y1;
xCir2 = r2 * cos(th) + x2; yCir2 = r2 * sin(th) + y2;

figure; hold on; grid on;
plot(xCir1, yCir1, 'color', 'b')
plot(xCir2, yCir2, 'color', 'r')

intercept_1 = [1.396562 1.752865]; % from code output above
intercept_2 = [-2.462136 4.968446]; 

scatter(intercept_1(1), intercept_1(2), 50, '+', 'MarkerEdgeColor', 'k');
scatter(intercept_2(1), intercept_2(2), 50, '+', 'MarkerEdgeColor', 'k');

xlim([-30, 20])
ylim([-30, 20])
