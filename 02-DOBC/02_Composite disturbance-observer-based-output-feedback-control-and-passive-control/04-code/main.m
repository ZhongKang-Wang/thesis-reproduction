clc; clear; close all;
%rng(42);
%% 仿真参数设置
h = 0.002; % 采样周期0.002s
T = 30; % 仿真总时长30s
t = 0: h: T - h; % 时间向量
N = length(t); % 时间向量长度

%% model 1
A1 = [2.2 -0.3; 0.1 -5];
F1 = [0.01; 0.01];
G1 = [-1.1; 0.1];
C1 = [1 0.6; 0.5 0.1];
D1 = [0.2 0.1];
H1 = [0.2; 0.1];
W1 = [0 1; -1 0];
V1 = [0.5 0];
M1 = [0.2; 0.4];

%% model 2
A2 = [1.2 -0.3; 2.1 -3];
F2 = [0.05; 0.01];
G2 = [-0.8; 0.2];
C2 = [2 0.6; 1.5 2.1];
D2 = [0.15 0.15];
H2 = [2; 1];
W2 = [0 1; -1 0];
V2 = [0.1 0];
M2 = [0.2; 0.2];

%% DOB
A1_hat = [-9.1078 -2.7607; 0.6636 -6.0615];
B1_hat = [-21.2809; 36.6158];
C1_hat = [-2.3135 1.2147];
D1_hat = 9.4092;
L1 = [14.6649; 4.6497];

A2_hat = [-29.8322 -2.2106; -3.4608 -32.2830];
B2_hat = [-191.6202; -151.2263];
C2_hat = [0.4170 -1.0293];
D2_hat = 41.6411;
L2 = [0.0557; 0.0386];

%% the transition probability matrix
Q = [-1 1; 2 -2];
P = expm(Q*h); % 转移概率矩阵
modes = zeros(1, N); % 模态矩阵

%% initial state
x = zeros(2, N);
x(:, 1) = [1.8; -0.5];

y = zeros(1, N);
y(:, 1) = D1*x(:, 1);

x_hat = zeros(2, N);

ew = zeros(2, N);
ew(:, 1) = [0.2; -0.1]; % the estimation error for the states of the disturbance d1(t)

gamma = 1;
lambda = 1;
% U = [0 0; 0 1];
z = zeros(2, N);
z(:, 1) = C1*x(:, 1);

w1_hat = zeros(2, N);
d1_hat = zeros(1, N); % 估计的外生扰动

w1 = zeros(2, N);
w1(:, 1) = ew(:, 1);

d1 = zeros(1, N);
d1(:, 1) = V1*w1(:, 1);

v = zeros(2, N);

%% 离线生成干扰
delta_t = sin(t)./(1 + t.^2); % 扰动δ(t)
d2 = 1 ./ (5 + 10.*t); % 干扰d2(t)

f = zeros(1, N);

%% 计算控制输入
u = zeros(1, N);

cur_mode = 1; % 初始模式
modes(1) = cur_mode;

cur_mode = markov_jump(P, cur_mode);
modes(2) = cur_mode;
% 获取当前模式参数
if cur_mode == 1
    A = A1; F = F1; G = G1; C = C1; D = D1; H = H1;
    W = W1; V = V1; M = M1;
    A_hat = A1_hat; B_hat = B1_hat; C_hat = C1_hat; D_hat = D1_hat;
    L = L1;
else
    A = A2; F = F2; G = G2; C = C2; D = D2; H = H2;
    W = W2; V = V2; M = M2;
    A_hat = A2_hat; B_hat = B2_hat; C_hat = C2_hat; D_hat = D2_hat;
    L = L2;
end

u(:, 1) = C_hat * x_hat(:, 1) + D_hat*y(:, 1) - d1_hat(:, 1);

%% 主循环
for k = 2: N

    d1(:, k) = V*w1(:, k-1);
    w1(:, k) = w1(:, k-1) + (W*w1(:, k-1) + M*delta_t(:, k-1))*h;
    
    y(:, k) = D*x(:, k-1); % the output measurement
    z(:, k) = C*x(:, k-1);

    x(:, k) = x(:, k-1) + (A*x(:, k-1) + F*f(:, k-1) + G*(u(:, k-1) + d1(:, k-1)) + H*d2(:, k-1))*h; 
    
    x_hat(:, k) = x_hat(:, k-1) + h*(A_hat*x_hat(:, k-1) + B_hat*y(:, k-1));

    v(:, k) = v(:, k-1) + ((W + L*D*G*V)*(v(:, k-1) - L*y(:, k-1)) + L*D*G*u(:, k-1))*h;

    d1_hat(:, k) = V*w1_hat(:, k-1);
    w1_hat(:, k) = v(:, k) - L*y(:, k);

    ew(:, k) = w1(:, k) - w1_hat(:, k);
    
    f(:, k) = x(2, k)*sin(t(k));
    
    cur_mode = markov_jump(P, cur_mode);
    modes(k+1) = cur_mode;

    % 更新控制输入
    if cur_mode == 1
        A = A1; F = F1; G = G1; C = C1; D = D1; H = H1;
        W = W1; V = V1; M = M1;
        A_hat = A1_hat; B_hat = B1_hat; C_hat = C1_hat; D_hat = D1_hat;
        L = L1;
    else
        A = A2; F = F2; G = G2; C = C2; D = D2; H = H2;
        W = W2; V = V2; M = M2;
        A_hat = A2_hat; B_hat = B2_hat; C_hat = C2_hat; D_hat = D2_hat;
        L = L2;
    end

    % 更新控制输入
    u(:, k) = C_hat * x_hat(:, k) + D_hat*y(:, k) - d1_hat(:, k);
    
end
%% display result
% switching signal
figure(2);
stairs(t, modes(1, 2:end), 'b', 'LineWidth', 1);
axis([0 30 0.5 2.5]);
xlabel('t/sec');
title('Switching signal');
print(gcf, 'Figure2', '-depsc2');  

% states of the system with the DOB output feedback controller
figure(3);
plot(t, x(1, :), 'r', 'LineWidth', 1.5);
hold on;
plot(t, x(2, :), 'g.', 'lineWidth', 1.5);
axis([0 30 -1.5 2]);
title('States of the system');
xlabel('t/sec');
legend('x1(t)', 'x2(t)');
print(gcf, 'Figure3', '-depsc2');  

% States of the output feedback controller
figure(4);
plot(t, x_hat(1, :), 'r', 'LineWidth', 1.5);
hold on;
plot(t, x_hat(2, :), 'g.', 'lineWidth', 1.5);
%axis([0 30 -3 4]);
xlabel('t/sec');
title('States of the output feedback controller');
legend('x-hat1(t)', 'x-hat2(t)');
print(gcf, 'Figure4', '-depsc2');  

% Estimation error of the state of the distance d1(t)
figure(5);
plot(t, ew(1, :), 'r', 'LineWidth', 1.5);
hold on;
plot(t, ew(2, :), 'g.', 'lineWidth', 1.5);
%axis([0 30 -3 4]);
xlabel('t/sec');
title('Estimation error ew(t) of the states of the disturbance d1(t)');
legend('ew1(t)', 'ew2(t)');
print(gcf, 'Figure5', '-depsc2');  

% d1(t)
figure(6);
plot(t, d1(1, :), 'LineWidth', 1.5);
hold on;
plot(t, d1_hat(1, :), 'LineWidth', 1.5);
plot(t, d1(1, :) - d1_hat(:, 1), 'LineWidth', 1.5);
hold off;
%axis([0 30 0 2]);
title('d1, d1-hat, d1-d1hat');
legend('d1(t)', 'd1-hat(t)', 'd1(t)-d1hat(t)');
print(gcf, 'Figure6', '-depsc2');  

figure(7);
plot(t, u(1, :), 'LineWidth', 1.5);
%axis([0 30 -5 10])
title('Control input u(t)');
print(gcf, 'Figure7', '-depsc2');  

figure(8);
plot(t, z(1, :), 'LineWidth', 1.5);
hold on;
plot(t, z(2, :), 'LineWidth', 1.5);
hold off;
%axis([0 30 0 2]);
title('Controlled output with DOBC and passive control');
legend('z1(t)', 'z2(t)');
print(gcf, 'Figure8', '-depsc2');  

%figure(9);
% 说的是在没有DOBC的情况
% 也就是控制策略不包含d1_hat的估计
% 我发现如果缺少对d1的补偿，那么输出会有些毛刺，并不是一条光滑曲线。


