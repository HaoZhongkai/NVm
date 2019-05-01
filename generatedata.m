global gama0 Ta Tb
gama0 = 10.7083e-4;
Ta = 5e3;
Tb = 5e3;
%%%总共产生800组数据，以一个矩阵写入文件，结果以另一个矩阵写入文件
%%%矩阵每组数据有1000个数据点，转换为行向量输入
%%%输出矩阵每组20个，依次为10个核的A,B值
%%%周围有1到10个核的情况各80组
B = [400,505,450,600];
N = [24,32,56,96];
fp1 = fopen('datain.txt','a+');
fp2 = fopen('dataout.txt','a+');
for i1 = 1:10
    wh0 = 1e-3*(80*rand(5,i1)+10);
    th0 = pi*rand(5,i1);
    t = (0.01:0.01:10)';
    PX = zeros(80,1000);
    imax = 4;
    jmax = 4;
    kmax = 5;
    k2 = 1;
    for i2 = 1:imax
        for j = 1:jmax
            for k1 = 1:kmax
                A0 = wh0(k1,:).*cos(th0(k1,:));
                B0 = wh0(k1,:).*sin(th0(k1,:));
                Px = generatePx(B(i2),N(j),A0,B0);
                PX(k2,:) = Px;
                AB = zeros(1,20);
                for l = 1:length(A0)
                    AB(2*l-1) =A0(l);
                    AB(2*l) = B0(l);
                end
                fprintf(fp1,'%f ',Px);
                fprintf(fp1,'\n');
                fprintf(fp2,'%f ',AB);
                fprintf(fp2,'\n');
                k2 = k2+1;
            end
        end
    end
end
fclose(fp1);
fclose(fp2);


function Px = generatePx(B0,N,A,B)
    global gama0 Ta Tb;
    t = (0.01:0.01:10)';
    wl = 2*pi*gama0*B0;
    beta = t*wl;
    M = ones(length(t),1);
    ep = 0.01*randn(length(t),1);
    for i = 1:length(A)
        wb = sqrt((2*pi*A(i)+wl)^2+(2*pi*B(i))^2);
        mz = (2*pi*A(i)+wl)/wb;
        mx = 2*pi*B(i)/wb;
        alpha = t*wb;
        Cphi = cos(alpha).*cos(beta)-mz*sin(alpha).*sin(beta);
        phi = acos(Cphi);
        n01 = mx^2*(1-cos(alpha)).*(1-cos(beta))./(1+Cphi);
        M = M.*(1-n01.*sin(N*phi/2).^2);
    end
    Px = 1/2*M.*exp(-2*N*t/Ta)+1/3+1/6*exp(-2*N*t/Tb)+ep;
    for i = 1:length(t)
        if Px(i)>1
             Px(i) = Px(i)-2*ep(i);
        end
    end
    Px = Px';
end

