%寻找在y中第n个数据点(可能为一个array)处附近N0个点的局部峰
function [y0,n0] = local_peak(y,n,N0)
y0 = zeros(1,length(n));
n0 = zeros(1,length(n));
for i = 1:length(n)
    [y0(i),n0(i)] = min(y(max(1,n(i)-N0):min(length(y),n(i)+N0)));
    n0(i) = n(i)-N0+n0(i)-1;
end