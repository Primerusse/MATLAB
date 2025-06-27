clc;
clear;
close all;
format short

%% 参数设置
a = 2;
l = 5;
F = 50 * 10 ^ 3;
E1I1 = 24 * 10 ^ 6;
E2A2 = 6 * 10 ^ 7;

%% 解析解计算
analytical_solution = - (5 * E2A2 * F * a ^ 6) / (3 * E1I1) / (6 * E1I1  * l + 4 * E2A2 * a ^ 3);

%% 离散化参数
dx = 0.005;
% 微元长度
df = 0.5;
% 迭代增加轴力的步长
nx = a / dx;
% 长度为a的微元数

%% 初始化边界条件
F_AB_end = df;
% 杆AB末端剪力初始猜测值
F_CD_end = F;
% 杆CD末端集中力

%% 迭代模拟
for i = 1 : 10 ^ 10
    % 末端尽量大，使充分迭代
    [F_AB, M_AB, W_AB, F_CD, M_CD, W_CD, diff_w, diff_L, w_B] = sim(F_AB_end, F_CD_end, l, F, E1I1, E2A2, dx, nx);
    if mod(i, 100) == 0
        % 进度条
        clc;
        disp('程序运行进度:')
        a = w_B / analytical_solution * 100;
        s = ['[', repmat('*', 1, floor(a / 4)), repmat('-', 1, 25 - floor(a / 4)), ']'];
        fprintf('%s', s)
        fprintf('%2.2f%%', a)
    end
    if abs(diff_w - diff_L) <= 0.000001
        % diff_w与diff_L近似相等
        clc;
        fprintf('算法值为: %.4f\n', w_B);
        fprintf('解析值为: %.4f\n', analytical_solution);
        fprintf('算法值与解析值相对误差为: %.4f%%\n', abs((w_B - analytical_solution) / analytical_solution) * 100);
        disp('杆AB的挠度矩阵:')
        disp(W_AB)
        disp('杆CD的挠度矩阵:')
        disp(W_CD)
        break
    end
    F_AB_end = F_AB_end + df;
end

%% 定义模拟函数
function [F_AB, M_AB, W_AB, F_CD, M_CD, W_CD, diff_w, diff_L, w_B] = sim(F_AB_end, F_CD_end, l, F, E1I1, E2A2, dx, nx)
    % 初始化杆件状态
    % 杆AB的剪力、弯矩、挠度矩阵
    F_AB = zeros(1, nx);
    M_AB = zeros(1, nx);
    W_AB = zeros(1, nx);
    % 杆CD的剪力、弯矩、挠度矩阵
    F_CD = zeros(1, 2 * nx);
    M_CD = zeros(1, 2 * nx);
    W_CD = zeros(1, 2 * nx);

    % 初始化边界条件
    F_AB(end) = F_AB_end;
    F_CD(end) = F_CD_end;

    % 剪力的传递
    % 杆AB
    for i = nx-1 : -1 : 1
        F_AB(i) = F_AB(i + 1);
    end
    % 杆CD
    for i = 2 * nx - 1 : - 1 : nx + 1
        F_CD(i) = F_CD(i + 1);
    end
    for i = nx - 1 : - 1 : 1
        F_CD(i) = F - F_AB(i);
    end

    % 弯矩的传递
    % 杆AB
    for i = nx - 1 : - 1 : 1
        M_AB(i) = sum(F_AB(i + 1 : end) * dx);
    end
    % 杆CD
    for i = 2 * nx - 1 : - 1 : 1
        M_CD(i) = sum(F_CD(i + 1 : end) * dx);
    end

    % 挠度的传递
    % 挠度的二阶导数
    W2_AB = - M_AB / E1I1;
    W2_CD = - M_CD / E1I1;
    % 第一次积分，辛普森法数值积分
    slope_AB = simpson_integrate(W2_AB, dx);
    slope_CD = simpson_integrate(W2_CD, dx);
    % 第二次积分，辛普森法数值积分
    W_AB = simpson_integrate(slope_AB, dx);
    W_CD = simpson_integrate(slope_CD, dx);

    % 输出结果
    diff_w = W_AB(nx) - W_CD(nx);
    diff_L = F_AB(nx) * l / E2A2;
    w_B = W_AB(nx);
end

%% 辛普森法数值积分
function y = simpson_integrate(f, h)
    % 输入:
    %   f - 导数在等距点上的值向量 [f0, f1, f2, ..., fN]
    %   h - 采样步长
    % 输出:
    %   y - 原函数在对应点上的值向量 [y0, y1, y2, ..., yN]

    N = length(f) - 1;  % 采样点数-1
    y = zeros(1, N + 1);  % 存储原函数y的值

    % 初始条件
    y(1) = 0;

    % 逐点积分计算原函数y
    for k = 1 : N
        % 计算原函数y(k+1)
        if mod(k, 2) == 1
            % 奇步长: 梯形法则
            y(k + 1) = y(k) + (h / 2) * (f(k) + f(k + 1));
        else
            % 偶步长: 辛普森法则
            y(k + 1) = y(k - 1) + (h / 3) * (f(k - 1) + 4 * f(k) + f(k + 1));
        end
    end
end
