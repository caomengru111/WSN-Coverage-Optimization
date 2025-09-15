%% DE-PSO混合算法: WSNs节点覆盖率优化
clc; clear; close all;

%% 网络参数
L = 50;                % 区域边长
n = 35;                % 节点个数
R = 5;                 % 通信半径
data = 0.8;            % 离散粒度

%% 粒子群参数
maxgen = 500;          % 最大迭代次数（可根据需要增加）
sizepop = 30;          % 种群规模
Wmax = 0.9; Wmin = 0.4;% 惯性权重动态范围

c1 = 1.5;              % 自我认知参数
c2 = 2.0;              % 社会认知参数
Vmax = 6;              % 速度最大值
Vmin = -6;             % 速度最小值
popmax = L;            % 位置最大值
popmin = 0;            % 位置最小值

F = 0.5;               % DE缩放因子
CR = 0.9;              % DE交叉概率

%% 初始化种群
empty_pop.Position = [];
empty_pop.V = [];
pop = repmat(empty_pop, sizepop, 1);

fitness = zeros(sizepop, 1);

for i = 1:sizepop
    pop(i).Position = rand(n, 2) * L;   % 随机位置
    pop(i).V = Vmin + (Vmax - Vmin) * rand(n, 2); % 随机速度
    fitness(i) = fun(pop(i).Position(:,1), pop(i).Position(:,2), L, R, data);
end

[fitnessgbest, bestindex] = max(fitness);
gbest = pop(bestindex).Position;
pbest = pop;                       
fitnesspbest = fitness;           

disp('初始位置：');
disp(gbest);
disp(['初始覆盖率：', num2str(fitnessgbest)]);

%% 绘制初始分布
figure;
for i = 1:n
    axis([0 L 0 L]); axis square;
    hold on;
    sita = 0:pi/100:2*pi;
    fill(gbest(i,1)+R*cos(sita), gbest(i,2)+R*sin(sita), 'y');
end
plot(gbest(:,1), gbest(:,2), 'r*');
title('初始节点位置及覆盖区域');
legend('覆盖区域','WSNs节点');
xlabel('x'); ylabel('y');

%% 主循环
zz = zeros(maxgen, 1);

for gen = 1:maxgen
    W = Wmax - (Wmax - Wmin) * (gen / maxgen); % 惯性权重
    
    for j = 1:sizepop
        % ====== DE操作 ======
        idxs = randperm(sizepop, 3);
        while any(idxs == j)  % 保证不重复
            idxs = randperm(sizepop, 3);
        end
        X1 = pop(idxs(1)).Position;
        X2 = pop(idxs(2)).Position;
        X3 = pop(idxs(3)).Position;
        
        V_DE = X1 + F * (X2 - X3);  % 差分变异
        V_DE = max(min(V_DE, popmax), popmin);
        
        % 交叉
        U = pop(j).Position;
        jrand = randi(n*2); % 至少一维交叉
        for d = 1:2*n
            if rand < CR || d == jrand
                row = ceil(d/2);
                col = mod(d-1,2)+1;
                U(row, col) = V_DE(row, col);
            end
        end
        
        fitU = fun(U(:,1), U(:,2), L, R, data);
        if fitU > fitness(j)
            pop(j).Position = U;
            fitness(j) = fitU;
        end
        
        % ====== PSO操作 ======
        pop(j).V = W * pop(j).V + ...
                   c1 * rand * (pbest(j).Position - pop(j).Position) + ...
                   c2 * rand * (gbest - pop(j).Position);
               
        pop(j).V = max(min(pop(j).V, Vmax), Vmin); % 限制速度

        pop(j).Position = pop(j).Position + pop(j).V; % 更新位置
        
        pop(j).Position = max(min(pop(j).Position, popmax), popmin); % 边界
        
        % 更新适应度
        fitness(j) = fun(pop(j).Position(:,1), pop(j).Position(:,2), L, R, data);
        
        % 更新个体最优
        if fitness(j) > fitnesspbest(j)
            pbest(j).Position = pop(j).Position;
            fitnesspbest(j) = fitness(j);
        end
    end
    
    % 更新全局最优
    [fitnessgbest, bestindex] = max(fitnesspbest);
    gbest = pbest(bestindex).Position;
    zz(gen) = fitnessgbest;
   
end

%% 绘制结果
disp('最优位置：');
disp(gbest);
disp(['最优覆盖率：', num2str(zz(end))]);

figure; plot(zz,'r','LineWidth',2);
title('覆盖率收敛曲线'); xlabel('迭代次数'); ylabel('覆盖率');

figure;
for i = 1:n
    axis([0 L 0 L]); axis square;
    hold on;
    sita = 0:pi/100:2*pi;
    fill(gbest(i,1)+R*cos(sita), gbest(i,2)+R*sin(sita), 'g');
end
plot(gbest(:,1), gbest(:,2), 'r*');
title('PSO优化后节点位置及覆盖区域');
legend('覆盖区域','WSNs节点');
xlabel('x'); ylabel('y');

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
        D = sqrt((m-x(i)).^2+(n-y(i)).^2);   % 计算坐标点到圆心的距离
        [m0, n0] = find(D <= R);             % 检测出圆覆盖点的坐标
        Ind = (m0-1).*col+n0;                % 坐标与索引转化
        M(Ind) = 1;                          % 改变覆盖状态
    end
    scale = sum(M(1:end))/(row*col);         % 计算覆盖比例
    z = scale;

    
end