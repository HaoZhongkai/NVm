%对于已经分离出峰的离散信号进行峰分离
%称为离散梳函数
%与处理原信号的梳函数不同，该函数仅算出在该点(可以是一个向量的点)的值
function y = Discrete_comb(loc,T0,delta)
N = floor((loc(end)/T0+1)/2);%调制的数量
y = zeros(length(loc),1);
for i = 1:N
    y = y+exp(-((loc-(2*i-1)*T0)/delta).^2);%高斯线型
%     y = y+1./(((t-(2*i-1)*T0)/delta).^2+1);
%         if (loc(j)-(2*i-1)*T0>=-delta) && (loc(j)-(2*i+1)*T0)<=delta
%             y(j) = y(j)+1;
%         end
end
