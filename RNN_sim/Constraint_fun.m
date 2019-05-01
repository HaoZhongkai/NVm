%由约束条件形成的偏差
%num0为约束的个数
function biase = Constraint_fun(wl,wb,param,num0,~)
    %计算Jacobi矩阵
    J = zeros(num0,length(param));
%     J = zeros(2*num0,length(param));
    len = (length(param))/2;
%     len = (length(param)-2*num0)/2;
    for i = 1:num0
        J(i,i) = 4*pi*(wl+2*pi*param(i));
        J(i,i+len) = 4*pi*param(i+len);
%         J(i,2*len+i) = 1;
%         J(i,2*len+num0+i) = 1;
    end
    %计算约束函数
%     h = zeros(2*num0,1);
%     h(1:num0) = (wl+2*pi*param(1:num0)).^2+(2*pi*param(len+1:len+num0)).^2-...
%         wb.^2*(1+bound_limit)^2+param(2*len+1:2*len+num0);
%     h(num0+1:2*num0) = wb.^2*(1-bound_limit)^2-(2*pi*param(1:num0)...
%         +wl).^2-(2*pi*param(len+1:len+num0)).^2+param(2*len+num0+1:end);
%     biase  = -J'/(J*J')*h;
    h = (wl+2*pi*param(1:num0)).^2+(2*pi*param(len+1:len+num0)).^2-wb.^2;
    biase = -J'/(J*J')*h;
end