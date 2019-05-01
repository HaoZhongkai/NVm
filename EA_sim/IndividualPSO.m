classdef IndividualPSO<matlab.mixin.Copyable
    properties
        N_c
        fitness
        param  %����,N_c*2
        best_param  %��ǰ���ӵ����Ž�
        best_fitness%��ǰ�������Ž�
        v      %����,ά�Ⱥ�paramһ��
    end
    methods
        function self = IndividualPSO(N_c,wh_max,S,fitfun)
            self.N_c = N_c;
            self.param = zeros(N_c,2);
            wh = 1e-3*wh_max*rand(N_c,1);
            th = pi*rand(N_c,1);
            self.v = abs(1e-3*wh_max-max(wh))*rand(N_c,2);
            A = wh.*cos(th);
            B = wh.*sin(th);
            self.param = [A,B];
            self.best_param = self.param;       %���ó�ʼ����ֵ
            self.best_fitness = 0;
            self.fitness = self.fitnessfun(S,fitfun);
        end
        
         function Init_param(self,param)
            for i = 1:min(self.N_c,length(param))                
                self.param(i,1) = param(i,1);
                self.param(i,2) = param(i,2);
            end
         end
         function Renew_param(self,w,param_g,c1,c2,wh_max)    
             %����λ�����ٶ�(��������Ӧֵ)
             %wΪȨֵ,V_gΪȫ�������ٶ�,c1,c2�ֱ�Ϊ����������ļ�Ȩ�ٶ�
             self.v = self.v*w+c1*rand()*(param_g-self.param)...
                 +c2*rand()*(self.best_param-self.param);
             self.param = self.param + self.v;
             %��֤���в����ڽ���
             out_index = find(abs(self.param)>wh_max);
             self.param(out_index) = 0;
             self.v(abs(self.v)>wh_max) = 0;
         end
         
         function fitness = fitnessfun(self,S,fitfun)
            params = self.param;
            wl = S.wl;
            N = S.N;
            t = S.t;
            wh = zeros(1,length(params));
            th = zeros(1,length(params));
            for i = 1:length(params)
                %B<0��Ӧ��Ϊ0
                if params(i,2)<0
                    fitness = 0;
                    return
                end
                wh(i) = sqrt(params(i,1)^2+params(i,2)^2);
                th(i) = atan(params(i,2)/params(i,1));
                if params(i,1)<0
                    th(i) = th(i)+pi;
                end
            end
            S1 = Kernal(wh,th,wl,N,t);
            S1.get_Px();
            fitness = fitfun.Get_fit(S1,S); %������Ӧ��
            if fitness>self.best_fitness   %��������ʱ������Ӧ�ȹ�ʽ
                self.best_param = self.param;
                self.best_fitness = fitness;
            end
         end
         
         
    end
end