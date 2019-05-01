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
%%%载入待拟合数据
y = load('LMtest.mat');
y = y.Px;

%%%生成用来探测峰的权重矩阵的列向量
W = Weight_Matrix(y);

%%%设置初始值,y为待拟合的数据
f_Beta = Cal_Px(t,Beta);
S = Cal_Square(y,f_Beta,W);
S0 = 0;
%%%记录最小损失
Smin = S;
Beta_best = Beta;
f_best = f_Beta;
gama = 1;
%%%用S0来跟踪上一轮的拟合损失，S来表示这一轮的拟合损失

I = eye(length(Beta));
for i = 1:100
    J = Jacobi2(t,Beta);
    %%%生成阻尼系数所配的矩阵
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

    
    