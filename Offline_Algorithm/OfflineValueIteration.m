function [P, K1, P_values, K1_values, iter] = OfflineValueIteration(A, B, Q1, K1, R, gamma, maxIter)
% record the value of each iteration   
P_values = zeros(3, 3, maxIter);
P = zeros(3, 3);
K1_values = zeros(1, 3, maxIter);
% the initial value of K1
K1_values(:, :, 1) = K1;
for iter = 1:maxIter
    % policy iteration
    % P_opt = dlyap(sqrt(gamma)*(A-B*K1)', Q1+K1'*R*K1);
    P_opt = Q1 + K1'*R*K1 + gamma*(A - B*K1)'*P*(A - B*K1);
    % policy update
    K1 = (R+gamma*B'*P_opt*B) \ (gamma*B'*P_opt*A);
    K1_values(:,:,iter + 1) = K1; 
    % check for convergence
    if norm(P - P_opt, 'fro') < 1e-8
        break;
    else
    % update P
        P = P_opt;
        P_values(:,:,iter + 1) = P;
    end
end

