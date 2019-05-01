%%%�����������ݼ�10*5000
%%%��1��10���˵�����������,ÿ�ζ��������ķ�
%%%�Ὣ2��������(B,N)�����������ǰ����Ϊ����
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
        B0 = 500*rand(num0,1)+500;       %B��500��1000֮��
        N =  2*randi([20,150],num0,1);    %40��300֮���ż��
        wh = 1e-3*(rand(num0,i)*80+10); %wh��Χ
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
            temp_data = [B0(i1)*1e-3,N(i1)/300,S.Px'];%�ֱ�ΪB,N,S.Px
            data(i1+(i-1)*num0,:) = temp_data;
            iterate = iterate+1;
            if mod(iterate,10)==0
                disp(iterate);
            end
        end
    end
end
%%%��A,B�ڵڶ���ά�Ƚ���������
function AB = catAB(A,B)
    AB = zeros(size(A,1),2*size(A,2));
    for i = 1:size(A,2)
        AB(:,2*i-1) = A(:,i);
        AB(:,2*i) = B(:,i);
    end
end



