function y = Addnoise(y0,e)
%����Ϊ���Ϊ1����̬�ֲ������,����ֻ��������
ep = e*randn(size(y0));
y = y0;
for i = 1:length(y0)
    y(i) = y0(i)+ep(i);
    if y(i)>=1
        y(i) = y(i)-2*ep(i);
    end
end