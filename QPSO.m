
%% 清空环境变量
clc; clear;

%% 网络参数
L = 50;   % 区域边长
n = 35;   % 节点个数
R = 5;    % 通信半径
data = 0.8; % 离散粒度

%% QPSO参数
maxgen = 500;          
sizepop = 30;           
beta = 1.5;             % QPSO收缩扩展因子
popmax = L;             
popmin = 0;

%% 初始化种群
gbest = zeros(n, 2);   
empty_pop.Position = [];
pop = repmat(empty_pop, sizepop, 1);
fitness = zeros(sizepop, 1);

for i=1:sizepop
    pop(i).Position = rand(n, 2).*L; 
    fitness(i) = fun(pop(i).Position(:,1), pop(i).Position(:,2), L, R, data);  
end

[fitnessgbest, bestindex] = max(fitness);
gbest = pop(bestindex).Position;    
pbest = pop;                        
fitnesspbest = fitness;             

%% 记录收敛曲线
zz = zeros(maxgen, 1);

%% QPSO迭代
for iter = 1:maxgen
    % 均值位置mbest
    mbest = zeros(n, 2);
    for k = 1:sizepop
        mbest = mbest + pbest(k).Position;
    end
    mbest = mbest / sizepop;

    for j = 1:sizepop
        P = (pbest(j).Position + gbest) / 2;

        % 更新每个节点坐标
        for dim = 1:2
            u = rand();
            if u == 0
                u = eps; % 防止log(0)
            end
            phi = rand();
            sign_val = 1;
            if phi < 0.5
                sign_val = -1;
            end
            pop(j).Position(:,dim) = P(:,dim) + ...
                sign_val * beta * abs(mbest(:,dim) - pop(j).Position(:,dim)) .* log(1/u);
        end

        % 边界处理
        pop(j).Position = max(pop(j).Position, popmin);
        pop(j).Position = min(pop(j).Position, popmax);

        % 适应度
        fitness(j) = fun(pop(j).Position(:,1), pop(j).Position(:,2), L, R, data);  
    end

    % 更新pbest和gbest
    for j = 1:sizepop
        if fitness(j) > fitnesspbest(j)
            pbest(j).Position = pop(j).Position;
            fitnesspbest(j) = fitness(j);
        end
        if fitness(j) > fitnessgbest
            gbest = pop(j).Position;
            fitnessgbest = fitness(j);
        end
    end
    zz(iter) = fitnessgbest;
end

%% 最终结果
disp('最优位置：');
disp(gbest);
disp(['最优覆盖率：', num2str(fitnessgbest)]);

figure; plot(zz, 'r', 'LineWidth', 2);
title('QPSO 收敛曲线');
xlabel('迭代次数'); ylabel('覆盖率');

% 绘制最终覆盖图
figure;
for i = 1:n
    axis([0 L 0 L]);
    sita = 0:pi/100:2*pi;
    hold on;
    p2 = fill(gbest(i,1)+R*cos(sita), gbest(i,2)+R*sin(sita), 'g');
end
p1 = plot(gbest(:,1), gbest(:,2), 'r*');
legend([p1, p2], {'WSNs节点', '覆盖区域'});
title 'QPSO优化后节点位置及覆盖区域';
xlabel('x'); ylabel('y');



function z = fun(x, y, L, R, data)
%% 覆盖率适应度函数
N = length(x);
[m, n] = meshgrid(0:data:L);
[row, col] = size(m);

% ⚡️ 这里必须初始化覆盖矩阵M
M = zeros(row, col);

for i = 1:N
    D = sqrt((m - x(i)).^2 + (n - y(i)).^2);
    [m0, n0] = find(D <= R);
    Ind = (m0 - 1) .* col + n0;
    M(Ind) = 1;
end
scale = sum(M(:)) / (row * col);
z = scale;
end
