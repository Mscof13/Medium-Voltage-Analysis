% Part 2: medium-voltage line driven by a sinusoidal voltage source (no load)
syms t

% Parameters:
R = 2.24;
L = 0.0169;
G = 0.0000191;
C = 0.000000382;
U_1max = 35000 * sqrt(2);
w = 2*pi*50;

% Initial conditions
u0 = 0;
du0 = 0;
tspan = linspace(0, 0.4, 100000);

% Analytical solution:
alpha = -(R*C + L*G)/(2*L*C);
beta = sqrt(4*(R*G + 1)*L*C - (R*C + L*G)^2)/(2*L*C);
A = (U_1max*(- C*L*w^2 + G*R + 1))/(C^2*L^2*w^4 + C^2*L^2*w^2 - ...
    2*C*G*L*R*w^2 - 2*C*L*w^2 + G^2*R^2 + 2*G*R + 1);
B = (C*L*U_1max*w)/(C^2*L^2*w^4 + C^2*L^2*w^2 - 2*C*G*L*R*w^2 - ...
    2*C*L*w^2 + G^2*R^2 + 2*G*R + 1);
U_2max = sqrt(A^2 + B^2);
phi = atan(B/A);
C_1 = u0 - U_2max;
C_2 = du0/B - (u0 - U_2max)*alpha/beta;
U_2 = exp(alpha*t)*(C_1*cos(beta*t) + C_2*sin(beta*t)) + ...
    U_2max*cos(w*t - phi);
U_2_numeric = matlabFunction(U_2);
U_2_values = U_2_numeric(tspan);

% Numerical solution:
x0 = [u0; du0];
p = [R, L, G, C, U_1max, w];
[t, y] = ode23(@(t, x) line_ode(t, x, p), tspan, x0);

% Plot:
figure;
hold on;
plot(tspan, U_2_values, 'r--', 'LineWidth', 3);
plot(tspan, y(:,1), 'b-', 'LineWidth', 1.5);
xlabel('$t [s]$', 'Interpreter', 'latex');
ylabel('$u_2(t) [V]$', 'Interpreter', 'latex');
title('Analytical and numerical solution', 'Interpreter', 'latex');
legend({'Analytical solution', 'Numerical solution'}, 'Interpreter', 'latex', ...
    'Location', 'best');
grid on;
hold off;

% Local function (must be at the end of a script file)
function dx = line_ode(t, x, p)
    R = p(1);
    L = p(2);
    G = p(3);
    C = p(4);
    U_1max = p(5);
    w = p(6);
    u_1 = U_1max * cos(w * t);
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = 1/(L*C)*(u_1 - (R*C + L*G)*x(2) - (R*G + 1)*x(1));
end