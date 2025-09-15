%% 清空环境变量
clc
clear

%% 网络参数
L = 50;                  % 区域边长
n = 35;                  % 节点个数
R = 5;                   % 通信半径
data = 0.8;                % 离散粒度

%% 蜻蜓算法参数
maxgen = 500;            % 迭代次数
sizepop = 30;            % 蜻蜓数量
gbest = zeros(n, 2);     % 最优解

%% 随机生成蜻蜓位置和适应度值
empty_pop.Position = [];
pop = repmat(empty_pop, sizepop, 1);
for i = 1:sizepop
    pop(i).Position = rand(n, 2) .* L;  % 初始蜻蜓位置
    fitness(i) = fun(pop(i).Position(:, 1), pop(i).Position(:, 2), L, R, data);  % 蜻蜓的适应度
end
[bestfitness, bestindex] = max(fitness);
gbest = pop(bestindex).Position;    % 群体最优极值
fitnessgbest = bestfitness;         % 种群最优适应度值

%% 初始结果显示
disp('初始位置：');
disp([num2str(gbest)]);
disp(['初始覆盖率：', num2str(fitnessgbest)]);

% 初始覆盖图
figure
for i = 1:n
    axis([0 L 0 L]);        % 限制坐标范围
    x = gbest(:, 1);
    y = gbest(:, 2);
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    p2 = fill(x(i) + R * cos(sita), y(i) + R * sin(sita), 'y');
end
p1 = plot(gbest(:, 1), gbest(:, 2), 'r*');
legend([p1, p2], {'WSNs节点', '覆盖区域'});
title '初始节点位置及覆盖区域';
xlabel('x');
ylabel('y');

%% 迭代寻优
for i = 1:maxgen
    for j = 1:sizepop
        % 随机选择邻居
        neighbors = randperm(sizepop, round(sizepop / 5));  % 随机选择20%个邻居
        bestNeighborFitness = -inf;
        bestNeighborPosition = pop(j).Position; % 初始化为当前蜻蜓位置
        
        % 选择适应度最高的邻居
        for k = 1:length(neighbors)
            neighborIndex = neighbors(k);
            if fitness(neighborIndex) > bestNeighborFitness
                bestNeighborFitness = fitness(neighborIndex);
                bestNeighborPosition = pop(neighborIndex).Position;
            end
        end
        
        % 更新位置
        pop(j).Position = pop(j).Position + rand * (bestNeighborPosition - pop(j).Position) + rand * (gbest - pop(j).Position);
        % 边界处理
        pop(j).Position = max(pop(j).Position, 0);
        pop(j).Position = min(pop(j).Position, L);
        
        % 更新适应度值
        fitness(j) = fun(pop(j).Position(:, 1), pop(j).Position(:, 2), L, R, data);
    end
    
    %% 群体极值更新
    [bestfitness, bestindex] = max(fitness);
    if bestfitness > fitnessgbest
        gbest = pop(bestindex).Position;
        fitnessgbest = bestfitness;
    end
    
    %% 每一代群体最优值存入 zz 数组
    zz(i) = fitnessgbest;
end

%% 结果显示
disp('最优位置：');
disp([num2str(gbest)]);
disp(['最优覆盖率：', num2str(zz(end))]);

%% 绘图
figure;
plot(zz, 'r', 'lineWidth', 2);          % 画出迭代图
title('算法过程', 'fontsize', 12);
xlabel('迭代次数', 'fontsize', 12);
ylabel('粒子覆盖率', 'fontsize', 12);

figure
for i = 1:n
    axis([0 L 0 L]);        % 限制坐标范围
    x = gbest(:, 1);
    y = gbest(:, 2);
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    p2 = fill(x(i) + R * cos(sita), y(i) + R * sin(sita), 'g');
end
p1 = plot(gbest(:, 1), gbest(:, 2), 'r*');
legend([p1, p2], {'WSNs节点', '覆盖区域'});
title 'Dragonfly优化后节点位置及覆盖区域';
xlabel('x');
ylabel('y');

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
    M = zeros(row, col);                % 初始化覆盖状态
    for i = 1:N
        D = sqrt((m - x(i)).^2 + (n - y(i)).^2);   % 计算坐标点到圆心的距离
        [m0, n0] = find(D <= R);             % 检测出圆覆盖点的坐标
        Ind = (m0 - 1) .* col + n0;                % 坐标与索引转化
        M(Ind) = 1;                          % 改变覆盖状态
    end
    scale = sum(M(1:end)) / (row * col);         % 计算覆盖比例
    z = scale;
end