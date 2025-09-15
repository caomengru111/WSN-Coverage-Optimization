%% 清空环境变量
clc
clear

%% 网络参数
L = 50;                  % 区域边长
n = 35;                  % 节点个数
R = 5;                   % 通信半径
data = 0.8;                % 离散粒度

%% WOA参数
maxgen = 500;            % 迭代次数
sizepop = 30;            % 鲸鱼个体数量
a = 2;                   % 控制系数a，线性递减

%% 参数初始化
gbest = zeros(n, 2);     % 全局最优位置

% 创建空的群体结构体
empty_pop.Position = [];
pop = repmat(empty_pop, sizepop, 1);

% 随机初始化群体位置及适应度值
fitness = zeros(sizepop, 1);
for i = 1:sizepop
    pop(i).Position = rand(n, 2) * L;  % 随机初始化位置
    fitness(i) = fun(pop(i).Position(:, 1), pop(i).Position(:, 2), L, R, data);  % 适应度
end

% 找到当前全局最优
[bestfitness, bestindex] = max(fitness);
gbest = pop(bestindex).Position;  % 全局最优位置
fitnessgbest = bestfitness;       % 全局最优适应度值

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
    p2 = fill(x(i) + R*cos(sita), y(i) + R*sin(sita), 'y');
end
p1 = plot(gbest(:, 1), gbest(:, 2), 'r*');
legend([p1, p2], {'WSNs节点', '覆盖区域'});
title('初始节点位置及覆盖区域');
xlabel('x');
ylabel('y');

%% 迭代寻优
for gen = 1:maxgen
    a = 2 - gen * (2 / maxgen);  % 线性递减系数a
    A = 2 * a * rand(sizepop, 1) - a;  % 计算A参数
    C = 2 * rand(sizepop, 1);          % 计算C参数
    
    for j = 1:sizepop
        r1 = rand();  % 随机数r1用于选择不同的更新策略
        D = abs(C(j) * gbest - pop(j).Position);  % 计算距离向量
        
        if r1 < 0.5  % 选择包围猎物或螺旋位置更新策略
            if abs(A(j)) < 1
                % 包围猎物更新位置
                pop(j).Position = gbest - A(j) * D;
            else
                % 随机搜索猎物
                X_rand = rand(n, 2) * L;  % 随机位置
                D_rand = abs(C(j) * X_rand - pop(j).Position);
                pop(j).Position = X_rand - A(j) * D_rand;
            end
        else
            % 螺旋位置更新，模拟鲸鱼螺旋式进攻猎物
            b = 1;        % 常数b定义螺旋形状
            l = -1 + rand();  % 在[-1, 1]范围内的随机数l
            dist_to_prey = abs(gbest - pop(j).Position);
            pop(j).Position = dist_to_prey * exp(b*l) .* cos(2*pi*l) + gbest;
        end
        
        % 边界处理
        pop(j).Position = max(pop(j).Position, 0);
        pop(j).Position = min(pop(j).Position, L);
        
        % 更新适应度值
        fitness(j) = fun(pop(j).Position(:, 1), pop(j).Position(:, 2), L, R, data);
        
        % 更新全局最优解
        if fitness(j) > fitnessgbest
            gbest = pop(j).Position;
            fitnessgbest = fitness(j);
        end
    end
    
    % 每一代群体最优值存入zz数组
    zz(gen) = fitnessgbest;
end

%% 结果显示
disp('最优位置：');
disp([num2str(gbest)]);
disp(['最优覆盖率：', num2str(zz(end))]);

%% 绘图
figure;
plot(zz, 'r', 'lineWidth', 2);  % 画出迭代图
title('WOA算法过程', 'fontsize', 12);
xlabel('迭代次数', 'fontsize', 12);
ylabel('覆盖率', 'fontsize', 12);

figure
for i = 1:n
    axis([0 L 0 L]);        % 限制坐标范围
    x = gbest(:, 1);
    y = gbest(:, 2);
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    p2 = fill(x(i) + R*cos(sita), y(i) + R*sin(sita), 'g');
end
p1 = plot(gbest(:, 1), gbest(:, 2), 'r*');
legend([p1, p2], {'WSNs节点', '覆盖区域'});
title('WOA优化后节点位置及覆盖区域');
xlabel('x');
ylabel('y');

%% 适应度函数：WSNs的覆盖率
function z = fun(x, y, L, R, data)
N = length(x);                      % 节点总个数
[m, n] = meshgrid(0:data:L);        % 离散化区域内的点
[row, col] = size(m);
M = zeros(row * col, 1);            % 初始化覆盖状态
for i = 1:N
    D = sqrt((m-x(i)).^2+(n-y(i)).^2);   % 计算坐标点到圆心的距离
    [m0, n0] = find(D <= R);             % 检测出圆覆盖点的坐标
    Ind = (m0-1).*col+n0;                % 坐标与索引转化
    M(Ind) = 1;                          % 改变覆盖状态
end
scale = sum(M(1:end))/(row*col);         % 计算覆盖比例
z = scale;
end