function [index,signal] = Period_search(t,y,N0,num)
%给定t,Px,N0为初始点的下标,num为搜索附近的点的数量,index为返回的下标值
len = length(t);
TN0 = N0;
n = floor((len/N0-1)/2)+1;
index = zeros(1,n);
signal = zeros(1,n);
i0 = 1;
while i0<=n
    N = floor((2*i0-1)*TN0);
    [y0,n0] = local_peak(y,N,num);
    index(i0) = n0;
    signal(i0) = y0;
    TN0 = (floor(N/(2*i0-1))+TN0)/2+0.5;
    i0 = i0+1;
end
    