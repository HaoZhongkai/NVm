%记得信号复原的步骤为:
%1.加中心信号 2.衰减调制 3.加噪声
%可以省略某些步骤，但是不要反
classdef Kernal<handle
    properties
        wh
        th
        wl
        N
        t
        Px
    end
    methods
        function self = Kernal(wh,th,wl,N,t)
                self.wh = wh;
                self.th = th;
                self.wl = wl;
                self.N = N;
                self.t = t;
        end
%         Px的get属性
       
%         function set.Px(self,Px)
%             self.Px = Px;
%         end
        %-----
      
        function A = getA(self)
            A = self.wh.*cos(self.th);
        end
        
        function B = getB(self)
            B = self.wh.*sin(self.th);
        end
        
        function mz = get_mz(self)
            A = 2*pi*self.wh.*cos(self.th);
            wb = self.get_wb();
            mz = (self.wl+A)./wb;
        end
        
        function mx = get_mx(self)
            B = 2*pi*self.wh.*sin(self.th);
            wb = self.get_wb();
            mx = B./wb;
        end
        
        function wb = get_wb(self)
            wb = sqrt((2*pi*self.wh.*cos(self.th)+self.wl).^2+(2*pi*self.wh...
                .*sin(self.th)).^2);
        end
        
        function M = get_M(self)
            wb = self.get_wb();
            mz = self.get_mz();
            mx = self.get_mx();
            beta = self.t*self.wl;
            M = ones(length(self.t),1);
            for i = 1:length(self.wh)
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
        
        %加衰减调制信号
        function Modulate(self,Ta,Tb)
            M = 2*self.Px-1;
            self.Px = 1/2*M.*exp(-2*self.N*self.t/Ta)+1/3+1/6*exp...
                (-2*self.N*self.t/Tb);
        end
        %加误差
        function Addnoise(self,e)
            ep = e*randn(size(self.Px));
            self.Px = min(self.Px+ep,1);
        end
        %加中心峰信号,n为中心核的个数,
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

                

                
                
            
        
        
        