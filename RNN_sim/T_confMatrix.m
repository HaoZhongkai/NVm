%求解对本问题的约束条件的
%param为A,B组
%num为约束个数
function T = T_confMatrix(wl,param,num)
%J为Jacobi矩阵
    J = zeros(num,length(param));
    len = (length(param)-num)/2;
    for i = 1:num
        J(i,i) = 4*pi*(wl+2*pi*param(i));
        J(i,i+len) = 2*pi*param(i+len);
        J(i,2*len+i) = 1;
    end
    T = eye(length(param))-J'/(J*J')*J;
end