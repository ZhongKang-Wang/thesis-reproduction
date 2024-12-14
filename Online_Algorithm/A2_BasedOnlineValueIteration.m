clc; clear; close all;
% original system
A = [-1 2; 2.2 1.7]; B = [2; 1.6]; C = [1 2]; D = 0;
% calculate the eigenvalues of the original system
poles_1 = eig(A);
% the original system is unstable.
% because the eigenvalues are located outside the unit-open disc. 
% initialize parameters
Q = 6; 
R = 1; 
gamma = 0.8; 
F = -1;
I = 1; 
C1 = [C, -I];
Q1 = C1' * Q * C1;
% construct the augmented system
T = zeros(3, 3); 
T(1:2, 1:2) = A; T(3, 3) = F;
B1 = [B; 0];
C_aug = [C, 0];
% number of iteration
iter = 200;
% online value iteration
[P_opt, K1_opt, P_values, iter] = OnlineValueIteration(T, B1, Q1, R, gamma, iter);
% draw
P_ij_values = zeros(iter, 9);
for i = 1:3
    for j = 1:3
        P_ij_values(:, (i - 1) * 3 + j) = squeeze(P_values(i,j,1:iter));
    end
end
x_values = 0:40;
figure;
plot(x_values, P_ij_values(1:41, :)', 'LineWidth', 2);
xlabel(' Convergence of the P matrix parameters to their optimal values');
ylabel('P matrix parameters');
xlim = 0:5:40;
ylim = -50:0:200;
legend('P11','P12','P13','P22','P23','P33');
grid on;
hold off;




