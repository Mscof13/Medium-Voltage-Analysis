% Part 1: medium-voltage line driven by a constant current source
syms u(t) t

% Parameters:
G = 0.0000191;
C = 0.000000382;
I_0 = 200;

% Initial condition
u0 = 0;
tspan = linspace(0, 1, 1000);

% Analytical solution:
U_2 = (u0 - I_0/G)*exp(-G/C * t) + 1/G * I_0;
U_2_numeric = matlabFunction(U_2);
U_2_values = U_2_numeric(tspan);

% Numerical solution:
dudt = @(t, u) -(G/C)*u + 1/C * I_0;
[t, y] = ode45(dudt, tspan, u0);

% Plot:
figure;
hold on;
plot(tspan, U_2_values, 'r--', 'LineWidth', 3);
plot(t, y, 'b-', 'LineWidth', 1.5);
xlabel('$t [s]$', 'Interpreter', 'latex');
ylabel('$u_2(t) [V]$', 'Interpreter', 'latex');
title('Analytical and numerical solution', 'Interpreter', 'latex');
legend({'Analytical solution', 'Numerical solution'}, 'Interpreter', 'latex', ...
    'Location', 'best');
grid on;
hold off;