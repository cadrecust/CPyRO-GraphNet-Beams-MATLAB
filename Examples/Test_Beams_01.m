clc 
clear all

fcu = 35;
b = 400;
h = 600;
L = 6000;
dL = 100;

W1 = 20;
W2 = 30;

Ac = b * h;
Ic = b * h^3 / 12;
Ec = (3.46 * sqrt(fcu) + 3.21) * 1e3;

EI = Ec * Ic;

% These are the mechanic elements with the FEM
[R, U, V, M] = MSFSFEMBeams(L, Ac, Ic, Ec, [0, L], [W1, W2], dL, [0, L], 0);

%% BCs of mechanical elements
MLeft = M(1, 1);
[Mmid, mp] = max(M(1, :));
Mright = M(2, end);

VLeft = V(1, 1);
Vright = V(1, end);

%% Data points to enforce BC
xMmid = (mp - 1) * dL;
x0L = [10, xMmid, L-10] ;
%% Create continuous analytical bending moment function
% For a beam with linearly varying distributed load: w(x) = W1 + (W2-W1)*x/L
% The shear force: V(x) = VLeft - ?w(x)dx = VLeft - W1*x - (W2-W1)*x^2/(2L)
% The bending moment: M(x) = MLeft + ?V(x)dx
W1 = W1;
W2 = W2;
% Define x vector for continuous evaluation
x_continuous = linspace(0, L, 1000);

% Analytical bending moment equation (corrected)
% M(x) = MLeft + VLeft*x - (W1*x^2)/2 - (W2-W1)*x^3/(6L)
M_analytical = MLeft + VLeft * x_continuous + (W1 * x_continuous.^2) / 2 + ((W2 - W1) * x_continuous.^3) / (6 * L);

% Evaluate at specific points for verification
M_at_x0L = MLeft + VLeft * x0L' + (W1 * x0L'.^2) / 2 + ((W2 - W1) * x0L'.^3) / (6 * L);

%% Plot comparison
figure(1)
plot(x_continuous, M_analytical, 'b-', 'LineWidth', 2)
hold on
plot([10, xMmid, L-10], M_at_x0L, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
xlabel('Position along beam (mm)', 'FontSize', 12)
ylabel('Bending Moment (N-mm)', 'FontSize', 12)
grid on
legend('Analytical M(x)', 'Enforced BC points', 'Location', 'best')
title('Bending Moment Distribution', 'FontSize', 14)
