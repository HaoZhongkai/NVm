%% ���ⲿ��������
%-------------��������(Ҫ��������Ԥ����t��Px)
B0 = 506;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
wh = [];
th = [];
N = 24;
S0 = Kernal(wh,th,wl,N,t);
S0.Px = Px;

%% �Ŵ��㷨�������
wh_max = 100;
N_c = 12;       %�˵ĸ���
N0 = 1000;       %��ʼȺ�����
% fit_param = sqrt(S0.Px'*S0.Px);
fit_param = 1;
fitfun = Fitnessfun('MSE',fit_param);
anneal_param = 0.003; %�˻����
Num_last = 40;    %ѡ���Ժ���������
Num_max = 400;      %��Ⱥ�������
mutate_risk = 0.2;  %����������
mutate_ratio = 0.1; %�������
inverse_risk = 0.2; %���嵹λ����

iterate_num = 100;
select_best_prob = 0.8;  %�������ѡ����ѡ�����ĸ���
group_num = 1;          %3����Ⱥ
Init_param = [0.0783,0.04;0.04,0.034;0.0323,0.031;-0.013,0.018;-0.0222,...
    0.029;0.0157,0.02];
%%   ��Ⱥ�ĳ�ʼ��
load Community      %���ļ�����

%  Community = Population(S0,N0,N_c,wh_max,fitfun);
%set initial value                                %�趨��ʼֵ
% Community.SetInit(Init_param);

%%   ��Ⱥ���ݻ�
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
%     %��ͼ
%     clf();
%     plot(t,S0.Px,'color','r');
%     plot(t,S_best.Px,'color','g');
end

    
%% ���ӻ���������
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
%% ������Ÿ��岢�۲�
figure
hold on;
plot(t,S0.Px,'color','r');
plot(t,S_best.Px,'color','g');
axis([t(1) t(end) 0 1]);
figure
plot(Community.fitness);

%% �洢��������
disp(A_best);
disp(B_best);
scatter(A_best,B_best);
grid on;
save('Community.mat','Community');