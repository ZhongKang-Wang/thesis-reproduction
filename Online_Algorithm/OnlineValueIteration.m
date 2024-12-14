function [P, K1, P_values, iter] = OnlineValueIteration(A, B, Q1, R, gamma, maxIter)
P = Q1; % initial P
P_values = zeros(3, 3, maxIter); % record P of each iteration;
P_values(:, :, 2) = P;
x0 = [5; -5; 5]; % initial state
K1 = [0.3 1.3 0.75]; % initial K1
N = 6; % step
x_dat = zeros(3, N); % record xk of each iteration
x_dat(:, 1) = x0;
x_next_dat = zeros(3, N);
rhok = zeros(N, 1); % response
Ak = zeros(N, N); % xk'*P*xk 
Bk = zeros(N, N); % xk+1'*P*xk+1
for iter = 2: maxIter
    for i = 1: N
        xk = x_dat(:, i);
        u = -K1*xk;
        rhok(i) = xk'*(Q1 + K1'*R*K1)*xk; % response
        Ak(i, :) = [xk(1)^2 2*xk(1)*xk(2) 2*xk(1)*xk(3) ...
                xk(2)^2  2*xk(2)*xk(3) xk(3)^2]; % xk'Pxk
    
        x_next = A*xk + B*u;
        x_next_dat(:, i) = x_next;
        Bk(i, :) = [x_next(1)^2 2*x_next(1)*x_next(2) 2*x_next(1)*x_next(3) ...
                x_next(2)^2  2*x_next(2)*x_next(3) x_next(3)^2];
        x_dat(:,i+1) = x_next; 
    end
    x_dat(:, 1) = x_next_dat(:, N);
    % value update
    [U, S, V] = svd(Ak - gamma*Bk); % 奇异值分解
    S_inv = pinv(S); % 计算伪逆
    x = V*S_inv*U'*rhok; % P = [x1 x2 x3; x2 x4 x5; x3 x5 x6]
    P_opt = [x(1) x(2) x(3); x(2) x(4) x(5); x(3) x(5) x(6)];
    % policy improvement
    K1 = (R+gamma*B'*P_opt*B)\(gamma*B'*P_opt*A);
    % check for convergence
    if norm(P - P_opt, 'fro') < 1e-5
            break;
    else
    % update P
        P = P_opt;
        P_values(:,:,iter + 1) = P; 
    end
end


