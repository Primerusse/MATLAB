clc;
clear;
close all;

%% 参数设置
Lx = 1;
% x方向长度
Ly = 1;
% y方向长度
T = 5;
% 总时间
c = 1;
% 波速

%% 离散化参数
nx = 100;
% x方向网格数
ny = 100;
% y方向网格数
nt = 2000;
% 时间步数
dx = Lx/nx;
% x方向步长
dy = Ly/ny;
% y方向步长
dt = T/nt;
% 时间步长

%% 稳定性条件检查
CFL = c*dt*sqrt(1/dx^2 + 1/dy^2);
if CFL >= 1
    warning('CFL条件不满足，可能导致数值不稳定！');
end

%% 初始化网格
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);
t = linspace(0, T, nt);

%% 初始化u(t,y,x)
u = zeros(nt, nx, ny);

%% 边界条件
for i = 1 : nx
    for j = 1 : ny
        x0 = x(i);
        y0 = y(j);
        u(1, i, j) = exp(- ((x0 - 0.5 * Lx) .^ 2 + (y0 - 0.5 * Ly) .^ 2) / 0.02);
    end
end
u(2, :, :) = u(1, :, :);
u(:, 1, :) = 0;
u(:, end, :) = 0;
u(:, :, 1) = 0;
u(:, :, end) = 0;

%% 差分方程迭代
for n = 3 : nt
    for i = 2 : nx - 1
        for j = 2 : ny - 1
            u(n, i, j) = (c * dt / dx) ^ 2 * (u(n - 1, i - 1, j) - 2 * u(n - 1, i, j) + u(n - 1, i + 1, j))...
                        + (c * dt / dy) ^ 2 * (u(n - 1, i, j - 1) - 2 * u(n - 1, i, j) + u(n - 1, i, j + 1))...
                        - u(n - 2, i, j) + 2 * u(n - 1, i, j);
        end
    end
end

%% 画图
figure('Position', [100, 100, 800, 600]);
% 创建一个新的图形窗口，并设置其在屏幕上的位置和大小，四个参数分别为：窗口左下角的 x 坐标，窗口左下角的 y 坐标，窗口宽度，窗口高度
for n = 1 : nt
    [X, Y] = meshgrid(x, y);
    u0 = [];
    for y0 = 1 : ny
        u0 = [u0; u(n, :, y0)];
    end
    %u0是将1*nx*ny的三维数组二维化的nx*ny数组
    surf(X, Y, u0, 'EdgeColor', 'none');
    % 'EdgeColor', 'none'：设置曲面网格线的颜色为无（即不显示网格线），使图像更平滑
    colormap(jet);
    % 设置图形的颜色映射方案为jet
    colorbar;
    % 在图形右侧添加一个颜色条，显示颜色与数值的对应关系
    zlim([-2 2])
    title(sprintf('二维波动方程 u(x,y,t) 的数值解 t = %.3f', n*dt))
    xlabel('位置 x');
    ylabel('位置 y');
    zlabel('位移 u(x,y,t)');
    view(30, 30);
    pause(0.001);
    drawnow limitrate;
    % 强制MATLAB立即更新图形
end