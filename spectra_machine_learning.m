Iz=0.5*[1,0;0,-1];
Ix=0.5*[0,1;1,0];
Iy=0.5*[0,-1i;1i,0];

% N=64;
t=0.01:0.01:9;
n0_x=1:length(t);n0_y=1:length(t);n0_z=1:length(t);
n1_x=1:length(t);n1_y=1:length(t);n1_z=1:length(t);
n0=zeros(length(t),1);
n1=zeros(length(t),1);
phi=1:length(t);
P_x=zeros(length(t),1);
data_list=zeros(length(t),9);
for B0=400:50:600
    for N=[8,16,24,32,40,48,56,64]


    %     B0=500;
        w0=2*pi*(B0*0.0001)*10.7083;
        H0=w0*Iz;

        for A=10:10:150
            for B=10:10:150
                A_hf=2*pi*A*10^(-3);
                B_hf=2*pi*B*10^(-3);
                w1=sqrt((w0+A_hf)^2+B_hf^2);
                H1=(A_hf+w0)*Iz+B_hf*Ix;

                for k=1:length(t)
                    V0=expm(-1i*t(k)*H0)*expm(-2i*t(k)*H1)*expm(-1i*t(k)*H0);
                    V1=expm(-1i*t(k)*H1)*expm(-2i*t(k)*H0)*expm(-1i*t(k)*H1);
                    U0=V0^(N/2);
                    U1=V1^(N/2); 

                    n0(k)=norm([trace(logm(V0)*1i*Ix),trace(logm(V0)*1i*Iy),trace(logm(V0)*1i*Iz)]);
                    n0_x(k)=trace(logm(V0)*1i*Ix)/n0(k);
                    n0_y(k)=trace(logm(V0)*1i*Iy)/n0(k);
                    n0_z(k)=trace(logm(V0)*1i*Iz)/n0(k);    
                    n1(k)=norm([trace(logm(V1)*1i*Ix),trace(logm(V1)*1i*Iy),trace(logm(V1)*1i*Iz)]);
                    n1_x(k)=trace(logm(V1)*1i*Ix)/n1(k);
                    n1_y(k)=trace(logm(V1)*1i*Iy)/n1(k);
                    n1_z(k)=trace(logm(V1)*1i*Iz)/n1(k);

                    % phi(k)=acos(trace(V0)); 
                    
                    % note:n0=n1=2*real_phi

                    P_x(k)=(1+(1-(1-n0_x(k)*n1_x(k)-n0_y(k)*n1_y(k)-n0_z(k)*n1_z(k))*(sin(N*n0(k)/2))^2))/2;
                end

        %         filename = 'N_32_B0_500Gs_A_10kHz_B_10_kHz.xlsx';
                filename = sprintf('%s%d%s%d%s%d%s%d%s','E:\ml_test_data2\N_',N,'_B0_',B0,'Gs_A_',A,'kHz_B_',B,'_kHz.xlsx');

                tau_name = 'tau';
                n0_x_name = 'n0_x';n0_y_name = 'n0_y';n0_z_name = 'n0_z';        
                n1_x_name = 'n1_x';n1_y_name = 'n1_y';n1_z_name = 'n1_z'; 
                phi_name = 'phi';P_x_name = 'Px';
                name_list = {tau_name,n0_x_name,n0_y_name,n0_z_name,n1_x_name,n1_y_name,n1_z_name,phi_name,P_x_name};

                sheet = 1;
                xlswrite(filename,name_list,sheet,'A1')       

                data_list(:,1)=t';
                data_list(:,2)=n0_x';
                data_list(:,3)=n0_y';
                data_list(:,4)=n0_z';
                data_list(:,5)=n1_x';
                data_list(:,6)=n1_y';
                data_list(:,7)=n1_z';
                data_list(:,8)=(2*n0)';
                data_list(:,9)=P_x';

                xlswrite(filename,data_list,sheet,'A2')
            end
        end

    end
end