%% 清空环境变量
clc; clear;

%% 网络参数
L = 50;                  % 区域边长
n = 35;                  % 节点个数
R = 5;                   % 通信半径
data = 0.8;                % 离散粒度

%% 遗传算法参数
maxgen = 500;            % 迭代次数
sizepop = 30;            % 种群规模
pcrossover = 0.8;        % 交叉概率
pmutation = 0.2;         % 变异概率
popmax = 50;             % 位置最大值
popmin = 0;              % 位置最小值

%% 参数初始化
empty_pop.Position = [];  
pop = repmat(empty_pop, sizepop, 1);  
fitness = zeros(sizepop, 1);  

% 随机生成种群位置和适应度
for i=1:sizepop
    pop(i).Position = rand(n, 2).*L;  % 初始种群位置
    fitness(i) = fun(pop(i).Position(:, 1), pop(i).Position(:, 2), L, R, data);  % 计算适应度
end

% 找到最优解
[bestfitness, bestindex] = max(fitness);
gbest = pop(bestindex).Position;    % 种群最优解

%% 初始结果显示
disp('初始位置：');
disp([num2str(gbest)]);
disp(['初始覆盖率：', num2str(bestfitness)]);

% 初始覆盖图
figure;
for i = 1:n
    axis([0 L 0 L]);        % 限制坐标范围
    x = gbest(:, 1); 
    y = gbest(:, 2); 
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on; 
    p2 = fill(x(i)+R*cos(sita), y(i)+R*sin(sita), 'y');
end
p1 = plot(gbest(:, 1), gbest(:, 2), 'r*');
legend([p1, p2], {'WSNs节点', '覆盖区域'});
title '初始节点位置及覆盖区域';
xlabel('x');
ylabel('y');

%% 遗传算法迭代寻优
zz = zeros(1, maxgen);  % 存储每代最优值
for i = 1:maxgen
    % 选择父代
    selectedPop = select(pop, fitness);
    
    % 交叉操作
    newPop = crossover(selectedPop, pcrossover, popmin, popmax);
    
    % 变异操作
    newPop = mutation(newPop, pmutation, popmin, popmax);
    
    % 适应度更新
    for j = 1:sizepop
        fitness(j) = fun(newPop(j).Position(:, 1), newPop(j).Position(:, 2), L, R, data);
    end
    
    % 找到当前代的最优解
    [currentBestFitness, bestindex] = max(fitness);
    if currentBestFitness > bestfitness
        bestfitness = currentBestFitness;
        gbest = newPop(bestindex).Position;
    end
    
    % 记录最优值
    zz(i) = bestfitness;
    pop = newPop;  % 更新种群
end

%% 结果显示
disp('最优位置：');
disp([num2str(gbest)]);
disp(['最优覆盖率：', num2str(zz(end))]);

% 画出迭代图
figure;
plot(zz, 'r', 'lineWidth', 2);
title('算法过程', 'fontsize', 12);
xlabel('迭代次数', 'fontsize', 12);
ylabel('覆盖率', 'fontsize', 12);

% 优化后节点位置及覆盖区域显示
figure;
for i = 1:n
    axis([0 L 0 L]);        % 限制坐标范围
    x = gbest(:, 1); 
    y = gbest(:, 2); 
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on; 
    p2 = fill(x(i)+R*cos(sita), y(i)+R*sin(sita), 'g');
end
p1 = plot(gbest(:, 1), gbest(:, 2), 'r*');
legend([p1, p2], {'WSNs节点', '覆盖区域'});
title 'GA优化后节点位置及覆盖区域';
xlabel('x');
ylabel('y');

%% 函数定义
% 适应度函数：WSNs的覆盖率
function z = fun(x, y, L, R, data)
    N = length(x);                      
    [m, n] = meshgrid(0:data:L);        
    [row, col] = size(m);               
    M = zeros(row, col);  % 初始化覆盖矩阵

    for i = 1:N
        D = sqrt((m-x(i)).^2+(n-y(i)).^2);   % 计算坐标点到圆心的距离
        [m0, n0] = find(D <= R);             % 检测出圆覆盖点的坐标
        Ind = (m0-1).*col+n0;                % 坐标与索引转化
        M(Ind) = 1;                          % 改变覆盖状态
    end
    scale = sum(M(1:end))/(row*col);         % 计算覆盖比例
    z = scale;
end

% 锦标赛选择
function selected = select(pop, fitness)
    sizepop = numel(pop);
    selected = pop;
    for i = 1:sizepop
        competitors = randperm(sizepop, 2);
        if fitness(competitors(1)) > fitness(competitors(2))
            selected(i) = pop(competitors(1));
        else
            selected(i) = pop(competitors(2));
        end
    end
end

% 交叉操作
function newPop = crossover(pop, pcrossover, popmin, popmax)
    sizepop = numel(pop);
    newPop = pop;
    for i = 1:2:sizepop-1
        if rand < pcrossover
            alpha = rand;
            newPop(i).Position = alpha*pop(i).Position + (1-alpha)*pop(i+1).Position;
            newPop(i+1).Position = alpha*pop(i+1).Position + (1-alpha)*pop(i).Position;
            newPop(i).Position = min(max(newPop(i).Position, popmin), popmax);
            newPop(i+1).Position = min(max(newPop(i+1).Position, popmin), popmax);
        end
    end
end

% 变异操作
function newPop = mutation(pop, pmutation, popmin, popmax)
    sizepop = numel(pop);
    newPop = pop;
    for i = 1:sizepop
        if rand < pmutation
            mutationIndex = randi([1, numel(newPop(i).Position)], 1);
            newPop(i).Position(mutationIndex) = newPop(i).Position(mutationIndex) + randn;
            newPop(i).Position = min(max(newPop(i).Position, popmin), popmax);
        end
    end
end