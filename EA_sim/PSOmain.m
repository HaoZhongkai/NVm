%% 从外部加载数据
%-------------产生数据(要求工作区中预先有t和Px)
B0 = 675;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
wh = [];
th = [];
N = 24;
S0 = Kernal(wh,th,wl,N,t);
S0.Px = Px;

%% 粒子群算法参数设计
wh_max = 100;
N_c = 10;
N0 = 100;
fit_param = 1;
fitfun = Fitnessfun('MIX',fit_param);
c1 = 2;
c2 = 2;
w_s = 1.0;
w_f = 0.5;
max_iter = 400;
group_num = 1;

%% 种群初始化
Community = PopulationPSO(S0,N0,N_c,wh_max,fitfun,c1,c2,w_s,w_f,max_iter);


%% 种群演化
running_fitness = [];
for loop = 1:max_iter
    Community.ReNew();
    [best_fit,best_index] = Community.evaluate_fitness();
    disp(best_fit);
    running_fitness = vertcat(running_fitness,best_fit);
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
figure
plot(Community.fitness);

%% 存储结果与输出
disp(A_best);
disp(B_best);
scatter(A_best,B_best);
grid on;
save('CommunityPSO.mat','Community');








