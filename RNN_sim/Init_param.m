%Ϊ�Ѿ��������wb����һ���ֵ,param_max����Ԫ�طֱ�Ϊ�����С���ֵ
%�����Ժ󻹿���ֱ��������B�ĳ�ֵ����ô�Ժ󽫿��ԸĽ��������
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