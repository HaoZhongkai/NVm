%%%产生大量数据集10*5000
%%%从1到10个核的情况随机产生,每次都加入中心峰
%%%会将2个超参数(B,N)放在数组的最前面作为输入
gama0 = 10.7083e-4;
num0 = 5000;
n_nuclear = 10;
t = (0.01:0.01:10)';
[test_data,test_target] = GetData(t,10,1000,gama0);


function [data,target] = GetData(t,n_nuclear,num0,gama0)
    data = zeros(num0*n_nuclear,2+length(t));
    target = zeros(num0*n_nuclear,20);w
    iterate = 1;
    for i = 1:n_nuclear
        B0 = 500*rand(num0,1)+500;       %B在500到1000之间
        N =  2*randi([20,150],num0,1);    %40到300之间的偶数
        wh = 1e-3*(rand(num0,i)*80+10); %wh范围
        th = pi*rand(num0,i);
        A = wh.*cos(th);
        B = wh.*sin(th);
        target(((i-1)*num0+1):i*num0,1:2*i) = catAB(A,B);
        for i1 = 1:num0
            wl = 2*pi*gama0*B0(i1);
            S = Kernal(wh(i1,:),th(i1,:),wl,N(i1),t);
            S.get_Px();
            S.AddCentralSignal(30,20e-3);
            S.Addnoise(0.03);
            temp_data = [B0(i1)*1e-3,N(i1)/300,S.Px'];%分别为B,N,S.Px
            data(i1+(i-1)*num0,:) = temp_data;
            iterate = iterate+1;
            if mod(iterate,10)==0
                disp(iterate);
            end
        end
    end
end
%%%将A,B在第二个维度交错连起来
function AB = catAB(A,B)
    AB = zeros(size(A,1),2*size(A,2));
    for i = 1:size(A,2)
        AB(:,2*i-1) = A(:,i);
        AB(:,2*i) = B(:,i);
    end
end



