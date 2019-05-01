%第三种峰分离方法
%先用findpeak将峰的位置选出,然后再对已经选出的峰进行变换
%基于peak_seperation1
%输入数据为已经去掉基线的数据
%核的基本参数
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
tmax = 10;
tstep = 0.001;
t = (tstep:tstep:tmax)';
%衰减周期,误差,脉冲数,中心信号参数,这里不用考虑衰减
Ta = 5e3;
Tb = 5e3;
e = 0.01;
N = 56; 
N_center = 60;
wh_center = 1e-3*20;
%核的参数
wh = 1e-3*[83.8,47,55,19,33,25.1];
th = pi/180*[21,30,54,133,132,51];
% wh = 1e-3*[50,40,76];
% th = pi/180*[32,46,23];
%计算核的参数
Px3 = 0;

S = Kernal(wh,th,wl,N,t);
S.get_Px();
S.AddCentralSignal(N_center,wh_center);
S.Addnoise(e);
Px = S.Px;

%---------------
%使用findpeaks寻峰
Px1 = 1-Px;
[height,loc,pk_w,pk_p] = findpeaks(Px1,t,'MinPeakProminence',0.1,...
   'Annotate','extents');
%MinPeakDistance',0.005


figure
hold on;
plot(t,Px1);
scatter(loc,height,'MarkerEdgeColor','r');

%---------
%使用变换函数
T0 = 0.5:0.0001:0.7;
delta = 0.01;
Px2 = zeros(length(loc),1);
for i = 1:length(T0)
    Px2(i) = height'*Discrete_comb(loc,T0(i),delta);
end

figure
plot(T0,Px2);


%--下一步--
%----------------
%现在从信号中选取出最强的一个,筛选出Px中对应的峰值
%为了robust保证,我们在周围稍微搜索一下
N_search = 5;
[~,max_index] = max(Px2);
T_signal = T0(max_index);
N_period = floor((t(end)/T_signal+1)/2);
Px_signal = zeros(N_period,1);
t_signal = zeros(N_period,1);
for i = 1:N_period
    T_i = (2*i-1)*T_signal;
    [Px_signal(i),index] = min(Px(floor(T_i/tstep)-N_search:floor...
        (T_i/tstep)+N_search));
    t_signal(i) = tstep*(index+floor(T_i/tstep)-N_search-1);
end
figure
hold on;
plot(t,Px);
scatter(t_signal,Px_signal,'MarkerEdgeColor','r');






