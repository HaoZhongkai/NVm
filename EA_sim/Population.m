classdef Population<handle
    properties 
        num  %��Ⱥ����
        individual
        S
        iter_num%����
        fitness  %fit array
        best_individual
        bestfit
        Fitnessfun
    end
    methods
        %��Ⱥ��ʼ��,N0Ϊ��ʼ��Ⱥ����,N_cΪÿ������Ĳ�������(�˸���)
        %SΪĿ���Ż�����
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
        
        %�趨��ʼֵ
        %���ѡȡһЩ�����帳��ʼֵ(С��1/8��Ⱥsize)
        function SetInit(self,param)
           n_init = randi(round(1*self.num/8));
           init_list = sort(randi(self.num,[1,n_init]));
           for i = init_list
              self.individual{i}.Init_param(param);
              self.individual{i}.Renew(self.S,self.Fitnessfun)
           end
           [self.bestfit,self.best_individual] = self.evaluate_fitness();
        end
        %�鿴ĳ�������Kernal
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
        %��Ⱥ�и���ѡ��
        %anneal_paramΪ�˻����,ʹ���˻���exp(-T/anneal_param)
        function Boltzmann_select(self,anneal_param,num_last)
            %������ѡ���ʱ��fit��һ��
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
                %��������
            end
            %ѡ�������+1
            self.iter_num = self.iter_num+1;
        end
        
        %���������ѡȡ
        %����sb���ʲ����Ĳ���1������del_list����һ��ѡ��
        %�����Ⱥ������һ��λnum_last
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
        %����
        %ʹ��һ�½����㷨
        %����ǰ���뾭��selection
        %��������Ŀǰ��size�����size֮�����ѡһ��ֵ��Ϊ����Ⱥ����
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
        %��һ�½���Ļ�����������Ƭ�ν���
        %���ǿ��ǵ������ֲ������е�һ���������Ը���
        %ratioΪsegment����ĸ���
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
        %������̣�������ת��λ�����,���һ������Ϊ�������ı���
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
        
        %������Ӧ��
        function [best_fitness,best_individual] = evaluate_fitness(self)
            self.fitness = zeros(self.num,1);
            for i = 1:self.num
                self.fitness(i) = self.individual{i}.fitness;
                %��total fitness��һ����Ӧ��
            end
            [best_fitness,best_individual] = max(self.fitness);
            self.best_individual = best_individual;
            self.bestfit = best_fitness;
        end
        
        %������ӻ����������A,B
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

%һ���ܹ���һ���ĸ�������1�ĺ���,size��һ����,��ʾ���ɵ�ά��
function series = random_pick(probability)
    switch 1
        case 1
            series = ceil(-rand+probability);
        case 2
            series = rand(size(1),size(2));
            series = arrayfun(@(x) ceil(-x+probability),series);
    end
end