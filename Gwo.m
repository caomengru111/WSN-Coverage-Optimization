%% 清空环境变量
clc
clear
axis equal
%% 网络参数
L = 50;                  % 区域边长
n = 25;                  % 节点个数
R = 5;                   % 通信半径
data = 0.8;                % 离散粒度
%% GWO算法参数
maxgen = 500;            % 迭代次数
sizepop = 30;            % 群体规模

%% 参数初始化
popmax = 50;             % 位置最大值
popmin = 0;              % 位置最小值

% α, β, δ狼初始化
alpha.Position = zeros(n, 2);
alpha.fitness = -inf;

beta.Position = zeros(n, 2);
beta.fitness = -inf;

delta.Position = zeros(n, 2);
delta.fitness = -inf;

% 初始化群体
empty_pop.Position = [];
pop = repmat(empty_pop, sizepop, 1);
fitness = zeros(sizepop, 1);

for i=1:sizepop
    pop(i).Position = rand(n, 2) * L;  % 初始种群位置
    fitness(i) = fun(pop(i).Position(:, 1), pop(i).Position(:, 2), L, R, data);  % 适应度
    % 更新α, β, δ狼
    if fitness(i) > alpha.fitness
        delta = beta;
        beta = alpha;
        alpha.Position = pop(i).Position;
        alpha.fitness = fitness(i);
    elseif fitness(i) > beta.fitness
        delta = beta;
        beta.Position = pop(i).Position;
        beta.fitness = fitness(i);
    elseif fitness(i) > delta.fitness
        delta.Position = pop(i).Position;
        delta.fitness = fitness(i);
    end
end

%% 初始结果显示
disp('初始位置：');
disp([num2str(alpha.Position)]);
disp(['初始覆盖率：', num2str(alpha.fitness)]);

% 初始覆盖图
figure
for i = 1:n
    axis([0 L 0 L]);        % 限制坐标范围
    x = alpha.Position(:, 1);
    y = alpha.Position(:, 2);
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    p2 = fill(x(i) + R * cos(sita), y(i) + R * sin(sita), 'y');
end
p1 = plot(alpha.Position(:, 1), alpha.Position(:, 2), 'r*');
legend([p1, p2], {'WSNs节点', '覆盖区域'});
title '初始节点位置及覆盖区域';
% 添加黑色方框
rectangle('Position', [0 0 L L], 'EdgeColor', 'k', 'LineWidth', 2);
% 移除坐标轴
axis off

%% 迭代寻优
for i = 1:maxgen
    a = 2 - i * (2 / maxgen);  % 线性下降权重a
    
    for j = 1:sizepop
        % 更新位置
        A1 = 2 * a * rand(n, 2) - a;
        C1 = 2 * rand(n, 2);
        D_alpha = abs(C1 .* alpha.Position - pop(j).Position);
        X1 = alpha.Position - A1 .* D_alpha;
        
        A2 = 2 * a * rand(n, 2) - a;
        C2 = 2 * rand(n, 2);
        D_beta = abs(C2 .* beta.Position - pop(j).Position);
        X2 = beta.Position - A2 .* D_beta;
        
        A3 = 2 * a * rand(n, 2) - a;
        C3 = 2 * rand(n, 2);
        D_delta = abs(C3 .* delta.Position - pop(j).Position);
        X3 = delta.Position - A3 .* D_delta;
        
        % 更新个体位置
        pop(j).Position = (X1 + X2 + X3) / 3;
        
        % 边界处理
        pop(j).Position = max(pop(j).Position, popmin);
        pop(j).Position = min(pop(j).Position, popmax);
        
        % 更新适应度
        fitness(j) = fun(pop(j).Position(:, 1), pop(j).Position(:, 2), L, R, data);
        
        % 更新α, β, δ狼
        if fitness(j) > alpha.fitness
            delta = beta;
            beta = alpha;
            alpha.Position = pop(j).Position;
            alpha.fitness = fitness(j);
        elseif fitness(j) > beta.fitness
            delta = beta;
            beta.Position = pop(j).Position;
            beta.fitness = fitness(j);
        elseif fitness(j) > delta.fitness
            delta.Position = pop(j).Position;
            delta.fitness = fitness(j);
        end
    end
    
    % 每代最优值存入zz数组
    zz(i) = alpha.fitness;
end

%% 结果显示
disp('最优位置：');
disp([num2str(alpha.Position)]);
disp(['最优覆盖率：', num2str(zz(end))]);

%% 绘图
figure;
plot(zz, 'r', 'LineWidth', 2);          %  画出迭代图
title('算法过程', 'FontSize', 12);
xlabel('迭代次数', 'FontSize', 12);
ylabel('覆盖率', 'FontSize', 12);

figure
% 添加网格
grid_size = L/5;  % 网格大小
for i = 0:5
    % 垂直线
    line([i*grid_size i*grid_size], [0 L], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
    % 水平线
    line([0 L], [i*grid_size i*grid_size], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold on
for i = 1:n
    axis([0 L 0 L]);        % 限制坐标范围
    x = alpha.Position(:, 1);
    y = alpha.Position(:, 2);
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    p2 = fill(x(i) + R * cos(sita), y(i) + R * sin(sita), 'g');
end
p1 = plot(alpha.Position(:, 1), alpha.Position(:, 2), 'r*');
legend([p1, p2], {'WSNs节点', '覆盖区域'});
title 'GWO优化后节点位置及覆盖区域';
% 添加黑色方框
rectangle('Position', [0 0 L L], 'EdgeColor', 'k', 'LineWidth', 2);
% 移除坐标轴
axis off
function z = fun(x, y, L, R, data)
%% 适应度函数：WSNs的覆盖率
% input：
% x        圆心横坐标
% y        圆心纵坐标
% L        区域边长
% R        通信半径
% data     离散粒度
% output:
% z        覆盖率

N = length(x);                      % 节点总个数
[m, n] = meshgrid(0:data:L);        % 离散化区域内的点
[row, col] = size(m);
M = zeros(row, col);                % 覆盖状态初始化
for i = 1:N
    D = sqrt((m - x(i)).^2 + (n - y(i)).^2);   % 计算坐标点到圆心的距离
    [m0, n0] = find(D <= R);             % 检测出圆覆盖点的坐标
    Ind = (m0 - 1) .* col + n0;          % 坐标与索引转化
    M(Ind) = 1;                          % 改变覆盖状态
end
scale = sum(M(:)) / (row * col);         % 计算覆盖比例
z = scale;
end