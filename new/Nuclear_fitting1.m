%对所有核进行拟合,版本1.0
%使用峰值点方法(利用峰值点筛选峰以及利用峰值点拟合核参数)逐个对核进行拟合
%考虑中心峰，不考虑误差和基线


%-------------产生数据
%--------------------------------------------------------------------------
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
tmax = 10;
tstep = 0.001;
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
% wh = 1e-3*[33];
% th = pi/180*[132];
A0 = wh.*cos(th);
B0 = wh.*sin(th);
%计算核数据
S = Kernal(wh,th,wl,N,t);
S.get_Px();
S.AddCentralSignal(N_center,wh_center);
% S.Addnoise(e);
Px = S.Px;
beta0 = S.get_wb./wl;
% figure
% hold on;
% plot(t,Px,'color','g');


%-----------------核参数拟合-----
%--------------------------------------------------------------------------
beta = [];
wh_best = [];
th_best = [];
A_best = [];
B_best = [];
B_max = 0.05;

[T_signal,t_signal,Px_signal] = Peak_seperate(t,Px);
% scatter(t_signal,Px_signal);
Px_signal = 1-Px_signal;
len = length(t_signal);
%最小二乘求Beta
beta1 = pi/(wl*T_signal)-1;
[B_fit,signal1] = Fit_B(B_max,wl,beta1,N,t_signal,Px_signal);
A_fit = (sqrt((beta1*wl)^2-(2*pi*B_fit)^2)-wl)/(2*pi);
wh_fit = sqrt(A_fit^2+B_fit^2);
th_fit = pi+atan(B_fit/A_fit);
wh_best = horzcat(wh_best,wh_fit);
th_best = horzcat(th_best,th_fit);
A_best = horzcat(A_best,A_fit);
B_best = horzcat(B_best,B_fit);
S_fit = Kernal(wh_fit,th_fit,wl,N,t);
S_fit.get_Px();

figure
plot(1:length(t_signal),t_signal);
figure
hold on;
plot(t,S_fit.Px,'color','r');
plot(t,Px,'color','g');
scatter(t_signal,1-Px_signal);
scatter(t_signal,1-signal1);
%----------输出结果---------
disp(B0(5));
disp(A0(5));
disp(B_best);
disp(A_best);







