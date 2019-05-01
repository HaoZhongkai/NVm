function W = Weight_Matrix(y)
%%%由于它是对角矩阵，所以只返回一个列向量
n = length(y);
W = ones(n,1);
dipdectlimit = 0.8;
for i = 1:n
    if y(i)<dipdectlimit
        W(i) = 5;
    end
end
