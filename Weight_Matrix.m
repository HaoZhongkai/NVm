function W = Weight_Matrix(y)
%%%�������ǶԽǾ�������ֻ����һ��������
n = length(y);
W = ones(n,1);
dipdectlimit = 0.8;
for i = 1:n
    if y(i)<dipdectlimit
        W(i) = 5;
    end
end
