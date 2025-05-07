function next_mode = markov_jump(P, current_mode)
% current_mode: 当前模态(上一时刻)
% P: 转移概率矩阵
% next_mode: 生成模态(当前时刻)

% 生成0-1的随机数
switch_num = rand;
    
% 根据当前模式进行判断
switch current_mode 
    case 1
        if switch_num < P(1, 1)
            next_mode = 1;
        else
            next_mode = 2;
        end
    case 2
        if switch_num < P(2, 1)
            next_mode = 1;
        else
            next_mode = 2;
        end
    otherwise
        error('当前模式必须是1、2或3。');
end
end    