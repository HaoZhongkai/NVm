%�����Ѿ�����������ɢ�źŽ��з����
%��Ϊ��ɢ�ắ��
%�봦��ԭ�źŵ��ắ����ͬ���ú���������ڸõ�(������һ�������ĵ�)��ֵ
function y = Discrete_comb(loc,T0,delta)
N = floor((loc(end)/T0+1)/2);%���Ƶ�����
y = zeros(length(loc),1);
for i = 1:N
    y = y+exp(-((loc-(2*i-1)*T0)/delta).^2);%��˹����
%     y = y+1./(((t-(2*i-1)*T0)/delta).^2+1);
%         if (loc(j)-(2*i-1)*T0>=-delta) && (loc(j)-(2*i+1)*T0)<=delta
%             y(j) = y(j)+1;
%         end
end
