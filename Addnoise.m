function y = Addnoise(y0,e)
%噪声为振幅为1的正态分布随机数,并且只加下噪声
ep = e*randn(size(y0));
y = y0;
for i = 1:length(y0)
    y(i) = y0(i)+ep(i);
    if y(i)>=1
        y(i) = y(i)-2*ep(i);
    end
end