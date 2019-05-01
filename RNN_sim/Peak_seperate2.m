%将peak_seperation3写成一个函数,所有的参数列表如下(大部分需要在函数内部改)
%1.搜索的周期上下限:T0,T1
%2.搜索周期的步长
%3.寻峰函数相关细节优化
%4.comb函数的种类及步长
%5.寻出周期后在附近搜索的点数
%6.这里假定了时间是等步长的
%num为选择出的核的个数
%本函数将对应的wb和初值B选出来
function [wb,B_fit] = Peak_seperate2(t,wl,N,Px,num,B_max)
%---------------
%使用findpeaks寻峰
tstep = t(2)-t(1);
Px1 = Px;
[height,loc,~,~] = findpeaks(Px1,t,'MinPeakProminence',0.1,...
   'Annotate','extents');
%MinPeakDistance',0.005
%---------
%使用变换函数
%在wl周围百分之30查找
step_size = 0.0001;
T0 = pi/(2*wl)*0.7:step_size:pi/(2*wl)*1.3;
delta = 0.01;
Px2 = zeros(length(loc),1);
for i = 1:length(T0)
    Px2(i) = height'*Discrete_comb(loc,T0(i),delta);
end
%----------------
%现在从信号中选取出最强的几个,筛选出对应的周期与wb
%peak_T是对应的周期,T_index是在峰peak_height数组中的下标
[peak_height,peak_T] = findpeaks(Px2,T0,'MinPeakProminence',1);
[~,T_index] = maxk(peak_height,min(num,length(peak_height)));
T0_signal = peak_T(T_index);
wb = (pi./T0_signal-wl)';

%计算合适的B
%N_search是在周围点寻找最大值的范围
N_search = 4;
B_fit = zeros(num,1);
N_period = floor((t(end)./T0_signal+1)/2);    %分别计算对应的周期数

for i = 1:length(T0_signal)
    %计算出对该wb对应的信号的位置
    P_signal = zeros(N_period(i),1);              %预分配P与t(整个t上的周期)
    t_signal = zeros(N_period(i),1);
    for j = 1:N_period(i)
        T_i = (2*i-1)*T0_signal(i);
        [P_signal(i),index] = max(Px1(floor(T_i/tstep)-N_search:floor...
            (T_i/tstep)+N_search));
        t_signal(i) = tstep*(index+floor(T_i/tstep)-N_search-1);
    end
    B_fit(i) = Fit_B(B_max,wl,wb(i)/wl,N,t_signal,P_signal);
end


