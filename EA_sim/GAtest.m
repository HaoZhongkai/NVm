%% 从外部加载数据
%-------------产生数据(要求工作区中预先有t和Px)
B0 = 506;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
wh = [];
th = [];
N = 24;
S0 = Kernal(wh,th,wl,N,t);
S0.Px = Px;

%% 遗传算法参数设计
wh_max = 100;
N_c = 12;       %核的个数
N0 = 1000;       %初始群体个数
% fit_param = sqrt(S0.Px'*S0.Px);
fit_param = 1;
fitfun = Fitnessfun('MSE',fit_param);
anneal_param = 0.003; %退火参数
Num_last = 40;    %选择以后的最大数量
Num_max = 400;      %种群最大数量
mutate_risk = 0.2;  %个体变异概率
mutate_ratio = 0.1; %变异比率
inverse_risk = 0.2; %个体倒位概率

iterate_num = 100;
select_best_prob = 0.8;  %随机排序选择中选择最大的概率
group_num = 1;          %3个种群
Init_param = [0.0783,0.04;0.04,0.034;0.0323,0.031;-0.013,0.018;-0.0222,...
    0.029;0.0157,0.02];
%%   种群的初始化
load Community      %从文件加载

%  Community = Population(S0,N0,N_c,wh_max,fitfun);
%set initial value                                %设定初始值
% Community.SetInit(Init_param);

%%   种群的演化
running_fitness = [];
for loop = 1:iterate_num
%       Community.Boltzmann_select(anneal_param,Num_last);
        Community.Stochastic_SortSelection(select_best_prob,Num_last);
%         Community.uni_Crossover(Num_max);
        Community.SSUcrossover(Num_max,0.4);
%         Community.Single_Param_Mutate(mutate_risk,mutate_ratio,36);
        Community.Mutate(mutate_risk,mutate_ratio);
        Community.Inverse(inverse_risk);
        [best_fit,best_index] = Community.evaluate_fitness();
        disp(Community.bestfit);
        running_fitness = vertcat(running_fitness,best_fit);
%     %画图
%     clf();
%     plot(t,S0.Px,'color','r');
%     plot(t,S_best.Px,'color','g');
end

    
%% 可视化进化过程
for i = 1:group_num
    [A,B] = Community.Get_AB();
end
S_best = Community.Get_item(best_index);
S_best.get_Px();
A_best = S_best.wh.*cos(S_best.th);
B_best = S_best.wh.*sin(S_best.th);
figure
hold on;
grid on;
scatter(A,B,'Marker','.');
% scatter(A0,B0,'Marker','o','SizeData',30);
scatter(A_best,B_best,'o','SizeData',20);
figure
plot(running_fitness);
%% 输出最优个体并观察
figure
hold on;
plot(t,S0.Px,'color','r');
plot(t,S_best.Px,'color','g');
axis([t(1) t(end) 0 1]);
figure
plot(Community.fitness);

%% 存储结果与输出
disp(A_best);
disp(B_best);
scatter(A_best,B_best);
grid on;
save('Community.mat','Community');