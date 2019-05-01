%遗传算法


%% 数据准备
%-------------产生数据
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
tmax = 10;
tstep = 0.01;
t = (tstep:tstep:tmax)';
%%衰减周期,误差,脉冲数,中心信号参数,这里不用考虑衰减
Ta = 5e3;
Tb = 5e3;
e = 0.01;
N = 56; 
N_center = 60;
wh_center = 1e-3*20;
%核的参数
wh = 1e-3*[83.8,47,55,19,33,25.1];
th = pi/180*[21,30,54,133,132,51];
A0 = wh.*cos(th);
B0 = wh.*sin(th);
S0 = Kernal(wh,th,wl,N,t);
S0.get_Px();
S0.AddCentralSignal(N_center,wh_center);

%% 遗传算法参数设计
wh_max = 40;
N_c = 20;       %核的个数
N0 = 1000;       %初始群体个数
% fit_param = sqrt(S0.Px'*S0.Px);
fit_param = 1;
fitfun = Fitnessfun('MSE',fit_param);
anneal_param = 0.003; %退火参数
Num_last = 40;    %选择以后的最大数量
Num_max = 400;      %种群最大数量
mutate_risk = 0.2;  %个体变异概率
mutate_ratio = 0.3; %变异比率
inverse_risk = 0.2; %个体倒位概率
iterate_num = 100;
select_best_prob = 0.8;  %随机排序选择中选择最大的概率
group_num = 1;          %3个种群
Init_param = [0.0783,0.04;0.0408,0.03;0.0323,0.028;-0.13,0.02;...
    -0.222,0.04;0.0156,0.016];
%%   种群的初始化
% load Community      %从文件加载

% Community = Population(S0,N0,N_c,wh_max,fitfun);
%多种群遗传算法
Community = cell(1,group_num);
for i = 1:group_num
    Community{i} = Population(S0,N0,N_c,wh_max,fitfun);%随机生成
    %set initial value                                %设定初始值
    Community{i}.SetInit(Init_param);
end

%%   种群的演化
running_fitness = [];
for loop = 1:iterate_num
    for j = 1:group_num
%       Community.Boltzmann_select(anneal_param,Num_last);
        Community{j}.Stochastic_SortSelection(select_best_prob,Num_last);
%       Community.uni_Crossover(Num_max);
        Community{j}.SSUcrossover(Num_max,0.4);
%         if mod(i,3)==1
%             GroupMigrate(Community{randi(3)},Community{randi(3)},10);
%         end
%       Community.Single_Param_Mutate(mutate_risk,mutate_ratio,36);
        Community{j}.Mutate(mutate_risk,mutate_ratio);
        Community{j}.Inverse(inverse_risk);
        [best_fit,best_index] = Community{j}.evaluate_fitness();
        S_best = Community{j}.Get_item(best_index);
        S_best.get_Px();
        disp(Community{j}.bestfit);
        running_fitness = vertcat(running_fitness,best_fit);
%     %画图
%     clf();
%     plot(t,S0.Px,'color','r');
%     plot(t,S_best.Px,'color','g');
    end
end
    
%% 可视化进化过程
A = [];
B = [];
for i = 1:group_num
    [A1,B1] = Community{1}.Get_AB();
    A = horzcat(A,A1);
    B = horzcat(B,B1);
end
A_best = S_best.wh.*sin(S_best.th);
B_best = S_best.wh.*cos(S_best.th);
figure
hold on;
scatter(A,B,'Marker','.');
scatter(A0,B0,'Marker','o','SizeData',30);
scatter(A_best,B_best,'o','SizeData',20);
figure
plot(running_fitness);
%% 输出最优个体并观察
figure
hold on;
plot(t,S0.Px,'color','r');
plot(t,S_best.Px,'color','g');
figure
plot(Community{1}.fitness);



save('Community.mat','Community');






