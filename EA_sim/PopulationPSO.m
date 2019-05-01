classdef PopulationPSO<handle
    properties
        num   %种群个数
        individual
        S
        iter_num  %当前迭代次数
        wh_max
        fitness   %适用函数
        best_individual %最优个体(最优个体指所有个体历史最优中最优的一个)
        bestfit
        Fitnessfun
        c1      %社会学习率
        c2      %个体学习率
        w_s     %初始比例
        w_f     %末态衰减比例
        max_iter
    end
    methods
        function self = PopulationPSO(S,N0,N_c,wh_max,fitfun,c1,c2,w_s,w_f...
                ,max_iter)
            self.individual = cell(1,N0);
            self.Fitnessfun = fitfun;
            self.wh_max = 1e-3*wh_max;
            for i = 1:N0
                self.individual{i} = IndividualPSO(N_c,wh_max,S,...
                    self.Fitnessfun);
            end
            self.num = N0;
            self.S = S;
            self.iter_num = 0;
            self.c1 = c1;
            self.c2 = c2;   %个体学习率与社会学习率
            self.w_s = w_s;
            self.w_f = w_f;
            self.max_iter = max_iter;
            [self.bestfit,self.best_individual] = self.evaluate_fitness();
        end
        
        %设定初始值
        %随机选取一些个个体赋初始值(小于1/8种群size)
        function SetInit(self,param)
           n_init = randi(round(1*self.num/8));
           init_list = sort(randi(self.num,[1,n_init]));
           for i = init_list
              self.individual{i}.Init_param(param);
              self.individual{i}.Renew(self.S,self.Fitnessfun)
           end
           [self.bestfit,self.best_individual] = self.evaluate_fitness();
        end
        
        function ReNew(self)
            w_k = self.w_f + (self.w_s-self.w_f)*self.iter_num/...
                self.max_iter;
            %所有个体调用参数更新公式
            param_g = self.individual{self.best_individual}.best_param;
            for i = 1:self.num   
                self.individual{i}.Renew_param(w_k,param_g,self.c1,...
                    self.c2,self.wh_max);
            end
        end
            
        
        function [best_fitness,best_individual] = evaluate_fitness(self)
            self.fitness = zeros(self.num,1);
            for i = 1:self.num
                self.fitness(i) = self.individual{i}.best_fitness;
                %用total fitness归一化适应度
            end
            [best_fitness,best_individual] = max(self.fitness);
            self.best_individual = best_individual;
            self.bestfit = best_fitness;
        end
        
        %方便可视化，提出所有A,B
        function [A,B] = Get_AB(self)
            A = zeros(1,self.num*self.individual{1}.N_c);
            B = zeros(1,self.num*self.individual{1}.N_c);
            for i = 1:self.num
                for j = 1:self.individual{1}.N_c
                    A((i-1)*self.individual{1}.N_c+j) = ...
                        self.individual{i}.param(j,1);
                    B((i-1)*self.individual{1}.N_c+j) = ...
                        self.individual{i}.param(j,2);
                end
            end
        end 
        
        %查看某个个体的Kernal
        function S = Get_item(self,index)
            len = length(self.individual{index}.param);
            wh = zeros(1,len);
            th = zeros(1,len);
            for i = 1:len
                param_i = self.individual{index}.param(i,:);
                wh(i) = sqrt(param_i(1)^2+param_i(2)^2);
                th(i) = atan(param_i(2)/param_i(1));
                if param_i(1)<0
                    th(i) = th(i)+pi;
                end
            end
            S = Kernal(wh,th,self.S.wl,self.S.N,self.S.t);
        end
    end
end