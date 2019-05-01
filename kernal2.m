classdef kernal
    properties 
        A
        B
    end
    methods
        function self = kernal(A,B)
            self.A = A;
            self.B = B;
        end
        function wb = get_wb(self)
            global wl;
            wb = sqrt((2*pi*self.A+wl)^2+(2*pi*self.B)^2);
        end
        function mx = get_mx(self)
            global wl;
            mx = (2*pi*self.B)/sqrt((2*pi*self.A+wl)^2+(2*pi*self.B)^2);
        end
        function mz = get_mz(self)
            global wl;
            mz = (2*pi*self.A+wl)/sqrt((2*pi*self.A+wl)^2+(2*pi*self.B)^2);
        end
        function a = get_alpha(self,tao)
            a = get_wb(self)*tao;
        end
        function p = phi(self,tao)
            global wl;
            p = acos(cos(self.get_alpha(tao)).*cos(wl*tao)...
                -self.get_mz()*sin(self.get_alpha(tao)).*...
                sin(wl*tao));
        end
        function M = get_M(self,tao)
            global N wl
            mx = self.get_mx();
            mz = self.get_mz();
            alpha = self.get_alpha(tao);
            beta = wl*tao;
            phi = self.phi(tao);
            M = 1-mx^2*(1-cos(alpha)).*(1-cos(beta))./(1+cos(phi)).*...
                (sin(N*phi/2).^2);
        end
        %%% 计算参数
        function m_der = get_m_der(self)
            mx = self.get_mx();
            mz = self.get_mz();
            wb = self.get_wb();
            m_der = [mx^2/wb,-mx*mz/wb;-mx*mz/wb,mz^2/wb];
        end
        %%%计算mz,mx分别对A,B的偏导数
        function alpha_der = get_alpha_der(self,tao)
            alpha_der = tao*[get_mz(self),get_mx(self)];
        end
        %%%计算alpha对A,B的导数，第一行是对A,第二行是对B
        function phi_der = get_phi_der(self,tao)
            global wl
            mz = get_mz(self);
            mx = get_mx(self);
            alpha = self.get_alpha(tao);
            beta = wl*tao;
            phi = self.phi(tao);
            wb = self.get_wb();
            phi_der = ((sin(alpha).*cos(beta)+mz*cos(alpha)...
                .*sin(beta)).*tao./sin(phi))*[mz,mx]-((mx*sin(alpha)...
                .*sin(beta))./(wb*sin(phi)))*[mx,-mz];
        end
        %%%计算phi对A,B的偏导数,第一行为对A,第二行为对B
        function M_der = get_M_der(self,tao)
            global N wl;
            mx = get_mx(self);
            alpha = self.get_alpha(tao);
            beta = wl*tao;
            phi = self.phi(tao);
            m_der = get_m_der(self);
            alpha_der = get_alpha_der(self,tao);
            phi_der = get_phi_der(self,tao);
            M1 = mx^2*sin(alpha).*(1-cos(beta)).*(sin(N*phi/2)...
                .^2)./(1+cos(phi));
            M2 = mx^2*(1-cos(alpha)).*(1-cos(beta)).*(sin(phi).*...
                (1-cos(N*phi))./(1+cos(phi))+N*sin(phi))./(2*(1+...
                cos(phi)));
            M_der =  2*mx*((1-cos(alpha)).*(1-cos(beta)).*(...
                sin(N*phi/2).^2)./(1+cos(phi)))*[m_der(2,1),m_der(2,2)]+...
                [alpha_der(:,1).*M1,alpha_der(:,2).*M1]+...
                [phi_der(:,1).*M2,phi_der(:,2).*M2];
        end
    end
end
    