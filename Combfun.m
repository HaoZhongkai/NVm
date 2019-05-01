%用周期性带有尖峰的函数排成一列调制信号
%T0为初始周期,delta大致为Gauss峰的宽度
%返回y为列向量
function y = Combfun(t,T0,delta)
N = floor((t(end)/T0+1)/2);%峰的数量
y = zeros(length(t),1);
for i = 1:N
    y = y+exp(-((t-(2*i-1)*T0)/delta).^2);%高斯线型
%     y = y+1./(((t-(2*i-1)*T0)/delta).^2+1);
end
