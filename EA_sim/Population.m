classdef Population<handle
    properties 
        num  %种群个数
        individual
        S
        iter_num%代数
        fitness  %fit array
        best_individual
        bestfit
        Fitnessfun
    end
    methods
        %种群初始化,N0为初始种群数量,N_c为每个个体的参数个数(核个数)
        %S为目标优化函数
        function self = Population(S,N0,N_c,wh_max,fitfun)
            self.individual = cell(1,N0);
            self.Fitnessfun = fitfun;
            for i = 1:N0
                self.individual{i} = Individual(N_c,wh_max,S,...
                    self.Fitnessfun);
            end
            self.num = N0;
            self.S = S;
            self.iter_num = 1;
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
        %种群中个体选择
        %anneal_param为退火参数,使用退火函数exp(-T/anneal_param)
        function Boltzmann_select(self,anneal_param,num_last)
            %仅需在选择的时候将fit归一化
            total_fit = 0;
            add_fit = 0;
            for i = 1:self.num
                add_fit = add_fit + self.individual{i}.fitness;
            end
            for i = 1:self.num
                total_fit = total_fit +exp(self.individual{i}.fitness...
                    /(add_fit*exp(-self.iter_num*anneal_param)));
            end
            %selection
            while(self.num>=num_last)
                i = randi(self.num);
                select_prob = exp(self.individual{i}.fitness/(add_fit*...
                    exp(-self.iter_num*anneal_param)))/...
                    total_fit*num_last;
                selection = random_pick(select_prob);
                if selection ~= 1
                    self.individual(i) = [];
                    self.num = self.num-1;
                end
                %致死个体
            end
            %选择完代数+1
            self.iter_num = self.iter_num+1;
        end
        
        %随机加排序选取
        %若以sb概率产生的不是1，则在del_list中挑一个选择
        %最后种群数量不一定位num_last
        function Stochastic_SortSelection(self,sb,num_last)
             [~,selection] = maxk(self.fitness,num_last);
             selection = selection';
             del_list = setdiff(1:self.num,selection);
             i = 1;
            while(i<=length(selection))
                if random_pick(sb) ~= 1
                    k = randi(length(del_list));
                    selection(i) = del_list(k);
                end
               i = i+1;
            end
            self.individual(setdiff(1:self.num,selection)) = [];
            self.num = length(self.individual);
            self.iter_num  = self.iter_num+1;
        end
        %交叉
        %使用一致交叉算法
        %交叉前必须经过selection
        %交叉后会在目前的size和最大size之间随机选一个值作为新种群数量
        function uni_Crossover(self,Max_num)
%             num_last = randi([self.num+10,Max_num]);
            num_last = Max_num;
            while(self.num <= num_last -2)
                i1 = randi([1,self.num]);
                i2 = randi([1,self.num]);
                chro1 = self.individual{i1}.chromosome;
                chro2 = self.individual{i2}.chromosome;
                [chro1,chro2] = uniform_crossover(chro1,chro2);
                new_indi1 = copy(self.individual{i1});
                new_indi2 = copy(self.individual{i2});
                new_indi1.chromosome = chro1;
                new_indi1.Renew(self.S,self.Fitnessfun);
                new_indi2.chromosome = chro2;
                new_indi2.Renew(self.S,self.Fitnessfun);
                self.individual = horzcat(self.individual,cell(1,2));
                self.individual{end-1} = new_indi1;
                self.individual{end} = new_indi2;
                self.num = self.num+2;
            end
        end
        %在一致交叉的基础上增加整片段交叉
        %这是考虑到将部分参数集中到一个个体明显改善
        %ratio为segment交叉的概率
        function SSUcrossover(self,Max_num,ratio)
            num_last = Max_num;
            while(self.num <= num_last -2)
                i1 = randi([1,self.num]);
                i2 = randi([1,self.num]);
                chro1 = self.individual{i1}.chromosome;
                chro2 = self.individual{i2}.chromosome;
                if rand>ratio
                    [chro1,chro2] = uniform_crossover(chro1,chro2);
                else
                    [chro1,chro2] = Segment_crossover(chro1,chro2);
                end
                new_indi1 = copy(self.individual{i1});
                new_indi2 = copy(self.individual{i2});
                new_indi1.chromosome = chro1;
                new_indi1.Renew(self.S,self.Fitnessfun);
                new_indi2.chromosome = chro2;
                new_indi2.Renew(self.S,self.Fitnessfun);
                self.individual = horzcat(self.individual,cell(1,2));
                self.individual{end-1} = new_indi1;
                self.individual{end} = new_indi2;
                self.num = self.num+2;
            end
        end
        
        function Single_Param_Mutate(self,mutate_risk,mutate_ratio,len)
            for i = 1:self.num
                mutate_prob = mutate_risk;
                if random_pick(mutate_prob)==1
                    mutate_seg = randi(self.individual{1}.N_c-1);
                    mutate_index = randi(self.num);
                    mutate_seg = (mutate_seg*len+1):(mutate_seg+1)*len;
                    chro_out = self.individual{mutate_index}.chromosome;
                    chro_out(mutate_seg) = mutate(chro_out(mutate_seg),...
                        mutate_ratio);
                    self.individual{mutate_index}.chromosome = chro_out;
                    self.individual{mutate_index}.Renew(self.S,self.Fitnessfun);
                end
            end
        end
        %变异过程，包括逆转和位点变异,最后一个参数为基因变异的比率
        function Mutate(self,mutate_risk,mutate_ratio)
            for i = 1:self.num
                mutate_prob = mutate_risk;
                if random_pick(mutate_prob)==1
                    mutate_index = randi(self.num);
                    chro_out = mutate(self.individual{mutate_index...
                        }.chromosome,mutate_ratio);
                    self.individual{mutate_index}.chromosome = chro_out;
                    self.individual{mutate_index}.Renew(self.S,self.Fitnessfun);
                end
            end
        end
        
        function Inverse(self,inverse_risk)
            for i = 1:self.num
                inverse_prob = inverse_risk;
                if random_pick(inverse_prob)==1
                    mutate_index = randi(self.num);
                    chro_out = inverse(self.individual{mutate_index...
                        }.chromosome);
                    self.individual{mutate_index}.chromosome = chro_out;
                    self.individual{mutate_index}.Renew(self.S,self.Fitnessfun);
                end
            end
        end
        
        %评估适应性
        function [best_fitness,best_individual] = evaluate_fitness(self)
            self.fitness = zeros(self.num,1);
            for i = 1:self.num
                self.fitness(i) = self.individual{i}.fitness;
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
    end
end

%一个能够以一定的概率生成1的函数,size是一个数,表示生成的维数
function series = random_pick(probability)
    switch 1
        case 1
            series = ceil(-rand+probability);
        case 2
            series = rand(size(1),size(2));
            series = arrayfun(@(x) ceil(-x+probability),series);
    end
end