function [T_signal,t_signal,Px_signal] = Peak_seperate_for_app(t,Px)
%---------------
%使用findpeaks寻峰
tstep = t(2)-t(1);
Px1 = 1-Px;
[height,loc,~,~] = findpeaks(Px1,t,'MinPeakProminence',0.1,...
   'Annotate','extents');
%MinPeakDistance',0.005


%---------
%使用变换函数
T0 = 0.5:0.0001:0.7;
delta = 0.01;
Px2 = zeros(length(loc),1);
for i = 1:length(T0)
    Px2(i) = height'*Discrete_comb(loc,T0(i),delta);
end

%----------------
%现在从信号中选取出最强的一个,筛选出Px中对应的峰值
%为了robust保证,我们在周围稍微搜索一下
N_search = 4;
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




