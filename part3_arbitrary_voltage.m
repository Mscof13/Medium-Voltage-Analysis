% Part 3: medium-voltage line driven by a piecewise-linear voltage source
% Analytical solution via the Laplace transform, numerical solution via ode23

% Initial conditions and source parameters
U_20 = 0;
U_20_prime = 0;
U_0 = 10;
t_1 = 1;
t_2 = 2;
t_3 = 4;

% Line parameters
R = 2.24;
L = 0.0169;
C = 0.000000382;
G = 0.0000191;

tspan = linspace(0, 10, 10000000);   % 10 million points; reduce for a faster run

% Analytical solution:
alpha = -(R*C + L*G)/(2*L*C);
beta = sqrt(4*(R*G + 1)*L*C - (R*C + L*G)^2)/(2*L*C);

a_1 = U_20/2;
b_1 = -(-U_20*alpha + U_20_prime)/(2*beta);

K_1 = -(U_0*t_1)/((t_2 - t_1)*(alpha^2 + beta^2));
a_2 = 1/2 * (U_0*t_1)/((t_2 - t_1)*(alpha^2 + beta^2));
b_2 = 1/2 * (U_0*t_1)/((t_2 - t_1)*(alpha^2 + beta^2))*alpha/beta;

K_2 = (U_0*t_2)/((t_2 - t_1)*(alpha^2 + beta^2));
a_3 = -1/2 * K_2;
b_3 = a_3*alpha/beta;

K_3 = -2*(U_0*alpha)/((t_2 - t_1)*(alpha^2 + beta^2)^2);
K_4 = K_3*(alpha^2 + beta^2)/(2*alpha);
a_4 = -1/2 * K_3;
b_4 = 1/2 * K_4/beta;

K_5 = -K_3;
K_6 = -K_4;
a_5 = -a_4;
b_5 = -b_4;

K_7 = 2*(U_0*alpha)/((t_3 - t_2)*(alpha^2 + beta^2)^2);
K_8 = K_7*(alpha^2 + beta^2)/(2*alpha);
a_6 = -1/2 * K_7;
b_6 = 1/2 * K_8/beta;

K_9 = -K_7;
K_10 = -K_8;
a_7 = -a_6;
b_7 = -b_6;

K_11 = -(U_0*t_3)/((t_3 - t_2)*(alpha^2 + beta^2));
a_8 = -1/2 * K_11;
b_8 = a_8*alpha/beta;

K_12 = (U_0*t_2)/((t_3 - t_2)*(alpha^2 + beta^2));
a_9 = -1/2 * K_12;
b_9 = a_9*alpha/beta;

U_2 = (2*a_1 .* exp(alpha.*tspan).*cos(beta.*tspan) - ...
    2*b_1.*exp(alpha.*tspan).*sin(beta.*tspan)) + ...
1/(L*C)*((K_2 + ...
2*a_3.*exp(alpha.*(tspan-t_1)).*cos(beta.*(tspan-t_1)) ...
- 2*b_3.*exp(alpha.*(tspan-t_1)).*sin(beta.*(tspan-t_1)) + ...
K_1 + 2*a_2.*exp(alpha.*(tspan-t_1)).*cos(beta.*(tspan-t_1)) - ...
2*b_2.*exp(alpha.*(tspan-t_1)).*sin(beta.*(tspan-t_1)) + ...
K_3 + K_4.*(tspan-t_1) + ...
2*a_4.*exp(alpha.*(tspan-t_1)).*cos(beta.*(tspan-t_1)) - ...
2*b_4.*exp(alpha.*(tspan-t_1)).*sin(beta.*(tspan-t_1))).*(tspan > t_1) + ...
(K_5 + K_6.*(tspan-t_2) + ...
2*a_5.*exp(alpha.*(tspan-t_2)).*cos(beta.*(tspan-t_2)) - ...
2*b_5.*exp(alpha.*(tspan-t_2)).*sin(beta.*(tspan-t_2)) + ...
K_7 + K_8.*(tspan-t_2) + ...
2*a_6.*exp(alpha.*(tspan-t_2)).*cos(beta.*(tspan-t_2)) - ...
2*b_6.*exp(alpha.*(tspan-t_2)).*sin(beta.*(tspan-t_2))).*(tspan > t_2) + ...
(K_9 + K_10.*(tspan-t_3) + ...
2*a_7.*exp(alpha.*(tspan-t_3)).*cos(beta.*(tspan-t_3)) - ...
2*b_7.*exp(alpha.*(tspan-t_3)).*sin(beta.*(tspan-t_3)) ...
+ K_11 + 2*a_8.*exp(alpha.*(tspan-t_3)).*cos(beta.*(tspan-t_3)) - ...
2*b_8.*exp(alpha.*(tspan-t_3)).*sin(beta.*(tspan-t_3)) ...
+ K_12 + 2*a_9.*exp(alpha.*(tspan-t_3)).*cos(beta.*(tspan-t_3)) - ...
2*b_9.*exp(alpha.*(tspan-t_3)).*sin(beta.*(tspan-t_3))).*(tspan > t_3));

% Numerical solution:
x0 = [U_20; U_20_prime];
p = [R, L, G, C, U_0, t_1, t_2, t_3];
[t, y] = ode23(@(t, x) line_ode(t, x, p), tspan, x0);

% Plot:
figure;
hold on;
plot(t, y(:, 1), 'b-', 'LineWidth', 1.5);
plot(tspan, U_2, 'r--', 'LineWidth', 2);
xlabel('$t [s]$', 'Interpreter', 'latex');
ylabel('$u_2(t) [V]$', 'Interpreter', 'latex');
title('Analytical and numerical solution', 'Interpreter', 'latex');
grid on;
hold off;

% Local function (must be at the end of a script file)
function dx = line_ode(t, x, p)
    R = p(1);
    L = p(2);
    G = p(3);
    C = p(4);
    U_0 = p(5);
    t_1 = p(6);
    t_2 = p(7);
    t_3 = p(8);

    % Source voltage: piecewise linear, drops from U_0 to 0 between t_1 and t_2,
    % rises back to U_0 between t_2 and t_3, and is zero otherwise
    excitation = ...
        (-U_0 / (t_2 - t_1) * t + U_0 / (t_2 - t_1) * t_2) * ((t >= t_1) & (t < t_2)) + ...
        (U_0 / (t_3 - t_2) * t - U_0 / (t_3 - t_2) * t_2) * ((t >= t_2) & (t < t_3));

    dx = zeros(2, 1);
    dx(1) = x(2);
    dx(2) = 1 / (L * C) * (excitation - (R * C + L * G) * x(2) - (R * G + 1) * x(1));
end