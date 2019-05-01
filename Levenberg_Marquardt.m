global wl N Ta Tb
B0 = 401;%magnetic field
gama0 = 10.7083e-4;%C13 parameter
wl = 2*pi*gama0*B0;
N = 56;
Tb = 5e3;
Ta = 5e3;
AB0 = [0.078234,0.030031,0.0407,0.0235];
Beta = 1.4e-2*randn(1,length(AB0)/2);
lambda = 0.02;
t = (0.01:0.01:10)';
%%%������������
y = load('LMtest.mat');
y = y.Px;

%%%��������̽����Ȩ�ؾ����������
W = Weight_Matrix(y);

%%%���ó�ʼֵ,yΪ����ϵ�����
f_Beta = Cal_Px(t,Beta);
S = Cal_Square(y,f_Beta,W);
S0 = 0;
%%%��¼��С��ʧ
Smin = S;
Beta_best = Beta;
f_best = f_Beta;
gama = 1;
%%%��S0��������һ�ֵ������ʧ��S����ʾ��һ�ֵ������ʧ

I = eye(length(Beta));
for i = 1:100
    J = Jacobi2(t,Beta);
    %%%��������ϵ������ľ���
    I = diag(diag(J'*diag(W)*J));
    delta = (J'*diag(W)*(y-f_Beta))\(J'*diag(W)*J+lambda*I);
    Beta = Beta+delta;
    f_Beta = Cal_Px(t,Beta);
    S0 = S;
    S = Cal_Square(y,f_Beta,W);
    gama = abs((S0-S)/S0);
    if S<S0
        lambda = lambda/1.2;
        if S<Smin
            Smin = S;
            Beta_best = Beta;
            f_best = f_Beta;
        end
    else
        lambda = lambda*1;
    end
end
figure
hold on;
plot(t,y,'color','r');
plot(t,f_best,'color',rand(1,3));

    
    