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
% start with a stablizing control policy K1
K1 = [0.3 1.3 0.75];
% number of iteration
iter = 200;
% offline policy iteration
[P_opt, K_opt, P_values, ~, iter] = OfflineValueIteration(T, B1, Q1, K1, R, gamma, iter);
% drawing
N = 50;  % step
x = zeros(3, N); 
y = zeros(1, N); 
u = zeros(1, N); 

% initial state
x(:,1) = [5; -5; 5];
y(1) = C_aug * x(:, 1);
for k = 1: N
    % calculate control input based on optimal control policy
    u(:, k) = -K_opt * x(:,k);
    x(:, k+1) = T * x(:, k) + B1 * u(:, k); % updata state
    y(k+1) = C_aug * x(:, k+1);
end
% 创建一个从0开始的 x 轴向量
x_values = 0: 40;
% plot the output changes
subplot(2, 1, 1);
plot(x_values, y(1, 1:41), 'g', 'LineWidth', 2);
hold on;
% plot the state r changes
subplot(2, 1, 1);
plot(x_values, x(3, 1:41), 'r', 'LineWidth', 2);
hold off;
xlabel('Time step');
ylabel('y and r');
title('Evaluation of the output and the reference trajectory');
legend('y', 'r');
% plot the input changes
subplot(2, 1, 2);
plot(x_values, u(:, 1:41), 'b', 'LineWidth', 2);
xlabel('Time step');
ylabel('Control input');
title('The control input during learning');
   

