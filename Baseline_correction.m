%本文件主要功能为基线的确定以及信号的去噪
%% 数据准备
%---------
%核的基本参数
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
tmax = 10;
tstep = 0.001;
t = (tstep:tstep:tmax)';
%衰减周期,误差,脉冲数,中心信号参数
Ta = 5e3;
Tb = 5e3;
e = 0.01;
N = 256; 
N_center = 50;
wh_center = 1e-3*20;
%核的参数
% wh = 1e-3*[83.8,47,55,19,33,25.1];
% th = pi/180*[21,30,54,133,132,51];
wh = 1e-3*[13];
th = pi/180*[75];
S = Kernal(wh,th,wl,N,t);
S.get_Px();
S.AddCentralSignal(N_center,wh_center);
S.Modulate(Ta,Tb);
S.Addnoise(e);
Px = S.Px;
%% 先使用原始数据进行一组基线拟合作为后面的初始值
fitfun = @(b,x)(1/2*exp(-b(1)*x)+1/6*exp(-b(2)*x)+1/3);
coef0 = nlinfit(t,Px,fitfun,[0,0]);

%% 使用cwt加sliding windows methods 去除尖峰
Px_diff = cwt(S.Px,0.0015*length(S.Px),'haar');
Px_diff([1:10,end-10:end]) = 0;
Px_diff = Px_diff';
%base_index为去除尖峰后数据的下标
[BaseIndex,PeakInfo] = SlidingWindows(Px_diff,8);
figure
plot(t,Px);
figure
plot(t,Px_diff);


%% 使用insymmetirc LSM方法求基线
%%%我们将直接求两个指数函数前的系数，再转过来求出Ta,Tb
weight = 0.01;       %默认的权重
rate = 0.00004;        %默认的梯度下降的速率
ep = 0.02;
%%%
t_data = t(BaseIndex);
Px_data = S.Px(BaseIndex);
coef = ASLSM_baseline(t,Px,fitfun,BaseIndex,ep,coef0,weight,rate);
Px_basefit = fitfun(coef,t);
%% 可视化与绘图
coef_0 = 2*N./[Ta,Tb];
% Baseline0 = fitfun(coef_0,t);
Baseline0 = 1/2*exp(-2*N*t/Ta)+1/6*exp(-2*N*t/Tb)+1/3;
figure
hold on;
plot(t,Px,'Marker','.');
plot(t,Px_basefit,'color','g','linewidth',2.5);
plot(t,Baseline0,'color','r');
% scatter(t(BaseIndex),S.Px(BaseIndex),'SizeData',1);
axis([0 tmax 0 1]);
T_fit = 2*N./coef;
disp(T_fit);






