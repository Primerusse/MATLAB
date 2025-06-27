clc;
clear; 
close all; 

%% 参数设置
L = 3;
% 空间范围 [0, L]
T = 6;
% 时间范围 [0, T]
c = 1;
% 波速

%% 离散化参数
nx = 100;
% 空间网格数
nt = 200;
% 时间网格数
dx = L / nx;
% 空间步长
dt = T / nt;
% 时间步长

%% CFL稳定性条件检查
CFL = c * dt / dx;
if CFL > 1
    warning('CFL条件不满足，可能导致数值不稳定！');
end

%% 初始化网格
x = linspace(0, L, nx);
t = linspace(0, T, nt);

%% 初始化u(t,x)，同一行表示同一时间，同一列表示同一空间
u = zeros(nt, nx);

%% 边界条件
u(1, :) = sin(2 * pi * x);
u(2, :) = u(1, :);
u(:, 1) = 0;
u(:, end) = 0;

%% 差分方程迭代
for n = 3 : nt
    for i = 2 : nx - 1
        u(n, i) = CFL ^ 2 * (u(n - 1, i - 1) - 2 * u(n - 1, i) + u(n - 1,i + 1)) - u(n - 2, i) + 2 * u(n - 1, i) + f(n, i);
    end
end

%% 绘制u(t,x)三维图形
[X, T] = meshgrid(x, t);
figure('Position', [100, 100, 800, 600]);
% 创建一个新的图形窗口，并设置其在屏幕上的位置和大小，四个参数分别为：窗口左下角的 x 坐标，窗口左下角的 y 坐标，窗口宽度，窗口高度
surf(X, T, u, 'EdgeColor', 'none');
% 'EdgeColor', 'none'：设置曲面网格线的颜色为无（即不显示网格线），使图像更平滑
colormap(jet);
% 设置图形的颜色映射方案为jet
colorbar;
% 在图形右侧添加一个颜色条，显示颜色与数值的对应关系
title('一维波动方程 u(x,t) 的数值解');
xlabel('位置 x');
ylabel('时间 t');
zlabel('位移 u(x,t)');
view(30, 30);  

%% 定义f(t, x)
function f=f(t, x)
    f = x;
end
