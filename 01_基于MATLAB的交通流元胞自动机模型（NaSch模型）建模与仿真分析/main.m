clc;
clear;
%输入参数
N = input('车辆数N:');
R = input('初始状态最右端车辆位置R：');
vmax = input('最大速度vmax:');
p = input('随机减速概率p:');
T = input('模拟迭代次数T:');

disp('初始状态：')
%随机生成1xN大小的递增的行向量xlt,要求所有元素大于等于1，小于等于R
xlt = sort(randperm(R, N))
%随机生成1xN大小的行向量vlt
vlt = randi([0, vmax], 1, N)

%进行T次迭代
for i = 1:T
    [xlt, vlt] = lteration(xlt, vlt, N, vmax, p);
end

%输出结果
disp('最终状态：')
xlt
vlt

%一次迭代函数定义
function [x, v] = lteration(xlt0, vlt0, n, vmax, p)
    x = xlt0;
    v = vlt0;
    x(n + 1) = inf;
    for i = 1:n
        v(i) = min(v(i) + 1, vmax);
        v(i) = min(x(i + 1) - x(i) - 1, v(i));
        % 以概率p执行
        if rand < p
            v(i) = max(v(i) - 1, 0);
        end
    end
    for i = 1:n
        x(i) = x(i) + v(i);
    end
    x(n + 1) = [];
end
