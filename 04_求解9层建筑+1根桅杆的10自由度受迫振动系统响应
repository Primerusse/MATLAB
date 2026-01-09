% 程序功能：求解9层建筑+1根桅杆的10自由度受迫振动系统响应
% 输出结果：各阶固有频率；总位移响应+各阶模态响应曲线
% 版本：v1.2

%% 环境初始化
clc; clear;
close all;
syms x(t) t;

%% 系统参数定义
% 建筑结构参数
n = 9;
% 楼层数量（总自由度=楼层数+桅杆=10）
m_build_sum = 1.395e9;
% 总楼重(kg)
m_build = m_build_sum / n;
% 每层建筑的质量(kg)
k_build_sum = 9.72e9;
% 楼层总刚度(N/m)
k_build = k_build_sum / n;
% 层间刚度(N/m)
c_build = 1.3e8;
% 层间阻尼系数(N·s/m)

% 桅杆结构参数（顶部附属结构）
m_mast = 5.76e6;
% 桅杆质量(kg)
k_mast = 8.54e9;
% 桅杆刚度(N/m)
c_mast = 6.6e5;
% 桅杆阻尼系数(N·s/m)

% 激励荷载参数
f = 1e7;
% 激励荷载幅值(N)
f_w = pi;
% 激励频率(rad/s)

% 瑞利阻尼参数（C=αM+βK）
alpha = 0.187;
% 质量比例阻尼系数
beta = 7.95e-5;
% 刚度比例阻尼系数

%% 系统矩阵构建
% 质量矩阵M
M = diag([m_build*ones(1,n), m_mast]);

% 刚度矩阵K
K = [2 * k_build, -k_build, zeros(1, n - 1)];
for i = 2: (n - 1)
    K = [K; zeros(1, i - 2), -k_build, 2 * k_build, -k_build, zeros(1, n - i)];
end
K = [K; zeros(1, n - 2), -k_build, k_build + k_mast, -k_mast; zeros(1, n - 1), -k_mast, k_mast];

% 阻尼矩阵C
C = alpha * M + beta * K;

% 荷载向量F
F = f * sin(f_w*t) * ones(n+1,1);

%% 模态分析
[phy, wn] = eig(K, M);
% phy：模态矩阵；wn：固有频率平方矩阵（对角元素为ω_i²）

% 振型归一化
% min(phy,[],1)：按列求最小值；./：矩阵按列相除
phy = phy ./ min(phy,[],1);
wn = sqrt(wn);

%% 模态解耦
MP = diag(phy' * M * phy);
KP = diag(phy' * K * phy);
CP = diag(phy' * C * phy);
FP = phy' * F;

%% 求解各阶模态响应（符号解）
sol_list = [];
for i = 1:n+1
    % 构建第i阶模态的单自由度振动方程：Mp(i) * x'' + Cp(i) * x' + Kp(i) * x = Fp(i)
    eqn = MP(i) * diff(x,t,2) + CP(i) * diff(x,t) + KP(i) * x == FP(i);

    % 初始条件：t=0时，位移=0，速度=0
    cond1 = x(0) == 0;
    cond2 = subs(diff(x,t), t, 0) == 0;

    % 求解符号微分方程并简化结果
    sol_i = simplify(dsolve(eqn, [cond1, cond2]));
    sol_list = [sol_list; sol_i];
end
sol_list = phy * sol_list;

%% 可视化
t_list = linspace(0, 100, 500);
cols = 4;
% 子图列数
rows = ceil((n + 1) / cols);
% 子图行数

figure('Name', '各阶模态位移曲线', 'Theme', "light", 'Position', [100, 100, 1200, 800]);
% 总标题
sgtitle(['各阶模态位移曲线（激励频率为', num2str(f_w), 'rad/s)'],'FontSize',14,'FontWeight','bold', 'Color','black');

for i = 1 : (n + 1)
    subplot(cols, rows, i);
    plot(t_list, vpa(subs(sol_list(i), t_list)), 'LineWidth',1.8, 'Color',[0.2,0.5,0.8]);
    wn_val = round(wn(i, i), 4);  % 提取第i阶固有频率（rad/s）
    title(['第', num2str(i), '阶（ω_n=',num2str(wn_val),'rad/s）'], 'FontSize',9,'FontWeight','bold', 'Color','black');
    grid on;
    grid minor;
end

%% 结果输出
% 输出各阶模态响应的符号解（便于理论分析）
disp('========================================');
disp('各阶固有频率（4位有效数字）：');
disp('========================================');
for i = 1:n+1
    disp(['第',num2str(i),'阶：',num2str(round(wn(i, i),4))]);
end
disp('========================================');

% 输出各阶模态响应的符号解（便于理论分析）
disp(' ');
disp('========================================');
disp('各阶响应（4位有效数字）：');
disp('========================================');
for i = 1:n+1
    disp(['第',num2str(i),'阶：',char(vpa(simplify(sol_list(i)),4))]);
end
disp('========================================');
