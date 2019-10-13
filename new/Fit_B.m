%ÄâºÏB
function [B_best,signal1] = Fit_B(B_max,wl,beta0,N,t_signal,Px_signal)
B_arr = linspace(0,B_max,100)';
Biase = zeros(length(B_arr),1);
for i = 1:length(B_arr)
    signal1 = Single_peakdepth(B_arr(i),wl,beta0,N,t_signal);
    Biase(i) = (Px_signal-signal1)'*(Px_signal-signal1);
end
[~,best_index] = min(Biase);
B_best = B_arr(best_index);
signal1 = Single_peakdepth(B_arr(best_index),wl,beta0,N,t_signal);
end