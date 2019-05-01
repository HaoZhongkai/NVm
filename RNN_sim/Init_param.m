%为已经分离出的wb整定一组初值,param_max两个元素分别为主峰和小峰的值
%鉴于以后还可以直接搜索出B的初值，那么以后将可以改进这个函数
function param = Init_param(wl,wb,B1,param_num,param_max)
    wh = param_max(2)*rand(param_num,1);
    th = pi*rand(param_num,1);
    A = wh.*cos(th);
    B = wh.*sin(th);
    B(1:length(wb)) = B1;
    A1 = (sqrt(wb.^2-(2*pi*B).^2)-wl)/(2*pi);
    A(1:length(wb)) = A1;
    param = [A;B];
end