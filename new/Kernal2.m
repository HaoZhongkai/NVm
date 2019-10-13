%�ǵ��źŸ�ԭ�Ĳ���Ϊ:
%1.�������ź� 2.˥������ 3.������
%����ʡ��ĳЩ���裬���ǲ�Ҫ��
%����A,B����
classdef Kernal2<handle
    properties
        A
        B
        wl
        N
        t
        Px
    end
    methods
        function self = Kernal2(A0,B0,wl,N,t)
            self.A = A0;
            self.B = B0;
            self.wl = wl;
            self.N = N;
            self.t = t;
        end
%         Px��get����
       
%         function set.Px(self,Px)
%             self.Px = Px;
%         end
        %-----
      
        function wh = get_wh(self)
            wh = sqrt((2*pi*self.A+self.wl).^2+(2*pi*self.B).^2);
        end
        
        function th = get_th(self)
            th = atan(self.B/(self.A*2*pi+self.wl));
        end
        
        function mz = get_mz(self)
            A0 = 2*pi*self.A;
            wb = self.get_wb();
            mz = (self.wl+A0)./wb;
        end
        
        function mx = get_mx(self)
            B0 = 2*pi*self.B;
            wb = self.get_wb();
            mx = B0./wb;
        end
        
        function wb = get_wb(self)
            wb = sqrt((2*pi*self.A+self.wl).^2+(2*pi*self.B).^2);
        end
        
        function M = get_M(self)
            wb = self.get_wb();
            mz = self.get_mz();
            mx = self.get_mx();
            beta = self.t*self.wl;
            M = ones(length(self.t),1);
            for i = 1:length(self.get_wh())
                alpha = self.t*wb(i);
                Cphi = cos(alpha).*cos(beta)-mz(i)*sin(alpha).*sin(beta);
                phi = acos(Cphi);
                n01 = mx(i)^2*(1-cos(alpha)).*(1-cos(beta))./(1+Cphi);
                M = M.*(1-n01.*sin(self.N*phi/2).^2);
            end
        end
        function get_Px(self)
            M = self.get_M();
            self.Px = (1+M)/2;
        end
        
        %��˥�������ź�
        function Modulate(self,Ta,Tb)
            M = 2*self.Px-1;
            self.Px = 1/2*M.*exp(-2*self.N*self.t/Ta)+1/3+1/6*exp...
                (-2*self.N*self.t/Tb);
        end
        %�����
        function Addnoise(self,e)
            ep = e*randn(size(self.Px));
            self.Px = min(self.Px+ep,1);
        end
        %�����ķ��ź�,nΪ���ĺ˵ĸ���,
        function AddCentralSignal(self,n,wh_max)
            M = self.get_M();
            wh_c = wh_max*rand(1,n);
            th_c = pi*rand(1,n);
            Px_center = Get_Px(self.t,wh_c,th_c,self.wl,self.N);
            M_center = 2*Px_center-1;
            M = M.*M_center;
            self.Px = (1+M)/2;
        end
    end
end

                

                
                
            
        
        
        