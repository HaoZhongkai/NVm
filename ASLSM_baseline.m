%% ʹ�÷ǶԳ���С���˷�������ߣ�Ŀ����ʹ�÷ǶԳ�Ȩ�ؽ����߸����ĵ������
% tΪ�����ʱ��,base_indexΪ���ڻ��ߵĵ���±�,epΪ���������
% ���ķ���Ϊ�ݶ��½�(����)
%weightΪ��
function coef = ASLSM_baseline(t,Px,fun,Base_index,ep,coef0,weight,rate)
    t_base = t(Base_index);
    Px_base = Px(Base_index);
    coef = coef0;
%     h = animatedline('color',rand(1,3));
    i = 0;
    pointnum = 1;
    pointnum0 = 0;
    while(abs(pointnum-pointnum0)/pointnum>0.000&&i<4000)
        Px_fit = fun(coef,t_base);
        pointnum0 = pointnum;   %�ϴε�ֵ
        pointnum = Lossfun(Px_base,Px_fit,weight,ep);
        grad = Grad_num(fun,Px_base,t_base,weight,coef,ep);
        coef = coef+rate*grad;
%         addpoints(h,i,pointnum);
%         pause(0.001);
        i = i+1;
    end
    %���Բ���
end
%�������ĵ�ĸ�������
function pointnum = Lossfun(Px_base,Px_fit,weight,ep)
    pointnum = sum((Px_base>=Px_fit).*exp((Px_fit-Px_base)/ep)+weight*(Px_base...
        <Px_fit).*exp((Px_base-Px_fit)/ep));
end
%�����ݶ�gradΪ������(Ϊ��ƥ��coef���е����)
function grad = Grad_num(fun,Px_base,t_base,weight,coef,ep)
    step = 1e-6;%�̶�����:�󵼵Ĳ���
    Px_fit0 = fun(coef,t_base);
    pointnum0 = Lossfun(Px_base,Px_fit0,weight,ep);
    coef1 = coef+[step,0];
    coef2 = coef+[0,step];
    Px_fit1 = fun(coef1,t_base);
    Px_fit2 = fun(coef2,t_base);
    pointnum1 = Lossfun(Px_base,Px_fit1,weight,ep);
    pointnum2 = Lossfun(Px_base,Px_fit2,weight,ep);
    grad = [pointnum1-pointnum0,pointnum2-pointnum0]/ep;
end



