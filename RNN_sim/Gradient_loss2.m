%���ɳڱ������ݶȼ���
%���������S��Ŀ��P0��L2�����¶Բ���A,B���ݶȺ���
%������һ����׼�źź�һ������õ�kernal
%���outputΪһ��������
function grad = Gradient_loss2(P0,param,t,wl,N,param_num)
    ep = 1e-8;
    tstep = t(2)-t(1);
    M0 = 1-2*P0;
    len = param_num;
    A = param(1:len);
    B = param(len+1:2*len);
    M1 = Get_M2(t,A,B,wl,N);
    grad = zeros(2*len,1);
    biase = (M1-M0)'*(M1-M0)*tstep;
    for i = 1:len
        A2 = A;
        B2 = B;
        A2(i) = A(i)+ep;
        B2(i) = B(i)+ep;
        M2 = Get_M2(t,A2,B,wl,N);
        M3 = Get_M2(t,A,B2,wl,N);
        grad(i) = ((M2-M0)'*(M2-M0)*tstep-biase)/(4*ep);
        grad(i+len) = ((M3-M0)'*(M3-M0)*tstep-biase)/(4*ep);
    end
end