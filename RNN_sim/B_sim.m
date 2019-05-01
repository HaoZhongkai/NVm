%使用自相关函数方法来做B的拟合
%% B的拟合
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
P0 = 1-S0.Px;
%% 使用自相关函数方法来求B的极小值
%寻峰算法
num0 = 6;
wb = Peak_seperate(t,P0,wl,num0);
B_max = 90e-3;
B_min = 10e-3;
B_search = linspace(B_min,B_max,300);
B1 = B_search;
Fit = zeros(length(B_search),1);
P0_norm = sqrt(P0'*P0);
for j = 1:length(B_search)
    B1 = B_search(j);
    A1 = (sqrt(wb(1).^2-(2*pi*B1).^2)-wl)/(2*pi);
    P1 = (1-Get_M2(t,A1,B1,wl,N))/2;
    Fit(j) = P1'*P0/(sqrt(P1'*P1)*P0_norm);
%     Fit(j) = (P1-P0)'*(P1-P0)*tstep;
end
plot(B_search,Fit);





