%%%�����������ݼ�10*5000
%%%��1��10���˵�����������,ÿ�ζ��������ķ�
%%%��������̶���B��N
%%%�˵ĸ�����ͬ��ʱ����Ϊһ��Ԫ�����鴢��,��ͬ������ǩ����ͬ
gama0 = 10.7083e-4;
num0 = 500;

n_nuclear = 10;
t = (0.01:0.01:10)';
%%������
B = 500;
N = 96;
str = ['hyper param B:',num2str(B),'N :',num2str(N)];
Data = GetData(t,10,num0,gama0,B,N);


function Data = GetData(t,n_nuclear,num0,gama0,B0,N)
    Data = cell(n_nuclear,2);   %��һάΪ�˵�����,�ڶ�άΪ��������֤��
    iterate = 1;
    for i = 1:n_nuclear
        data = zeros(num0,length(t));   %��ÿ���˵���������data
        wh = 1e-3*(rand(num0,i)*80+10); %wh��Χ
        th = pi*rand(num0,i);
        A = wh.*cos(th);
        B = wh.*sin(th);
        Data{i,2} = catAB(A,B);         %label
        for i1 = 1:num0
            wl = 2*pi*gama0*B0;
            S = Kernal(wh(i1,:),th(i1,:),wl,N,t);
            S.get_Px();
            S.AddCentralSignal(30,20e-3);
%             S.Addnoise(0.03);
            temp_data = S.Px';      %ֻ��Px
            data(i1,:) = temp_data;
            iterate = iterate+1;
            if mod(iterate,10)==0
                disp(iterate);
            end
        end
        Data{i,1} = data;
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



