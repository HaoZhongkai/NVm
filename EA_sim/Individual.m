%遗传算法
%将精确到4位小数的参数A,B乘10000取整再转化为二进制串
%最多考虑6位数字,使用gray编码储存
%每一代的基因存储在一个3维logical数组N_c*18*2中
%参数也存储在2维数组中，N_c*2
%初始设定N个核(前几个核通过寻峰算法只设一个参数,其他核随机化参数)
%gray码使用数组方便转为二进制串
%基因，染色体均为logical array
classdef Individual<matlab.mixin.Copyable
    properties
        gene   %基因(分开的数组)
        param        %参数,数组,N_c*2
        chromosome  %染色体,为将gene连接成一个字符串,为行向量
        N_c
        fitness
    end
    methods
        
        %初始化
        %S对应object function
        function self = Individual(N_c,wh_max,S,fitfun)
            self.N_c = N_c;
            self.param = zeros(N_c,2);
%             wh_pre = wh_max*rand(1,N_c);
            %核的初始值不断减小
            wh = 1e-3*wh_max*rand(1,N_c);
            th = pi*rand(1,N_c);
            A = wh.*cos(th);
            B = wh.*sin(th);
            for i=1:N_c
                self.param(i,:) = [A(i),B(i)];
            end
            self.param2gene();
            self.connact_gene();
            self.fitness = self.fitnessfun(S,fitfun);
        end
        
        function Init_param(self,param)
            for i = 1:min(self.N_c,length(param))                
                self.param(i,1) = param(i,1);
                self.param(i,2) = param(i,2);
            end
            self.param2gene();
            self.connact_gene();
        end
        
        %在染色体更新后对基因和参数进行更新
        %在产生产生给定条件的新个体时用到
         function Renew(self,S,fitfun)
            self.segment_chro();
            self.gene2param();
            self.fitness = self.fitnessfun(S,fitfun);
         end
         
         %将更新后染色体分割为基因
          function segment_chro(self)
              gene_length = length(self.chromosome)/(2*self.N_c);
             for i = 1:self.N_c
                 self.gene(i,:,1) = self.chromosome(...
                     2*(i-1)*gene_length+1:(2*i-1)*gene_length)';
                 self.gene(i,:,2) = self.chromosome(...
                     (2*i-1)*gene_length+1:2*i*gene_length)';
             end
         end
        
         function chro = connact_gene(self)
             chro = false(1,numel(self.gene));
             %每个核对应的基因长度，应为36
             len = length(chro)/self.N_c;
             for i = 1:self.N_c
                 chro((i-1)*len+1:i*len) = [self.gene(i,:,1)',...
                     self.gene(i,:,2)'];
             end
             self.chromosome = chro;
         end
         
        %输入param为一个二维数组(行)
        %gene为一个N_c*n*2的矩阵，代表A,B两个参数
        function gene = param2gene(self)
            gene = false(self.N_c,18,2);
            for i = 1:length(self.param)
                A = floor(self.param(i,1)*1e6);
                B = floor(self.param(i,2)*1e6);
                A_bi = de2bi(abs(A),18,'left-msb')';
                B_bi = de2bi(abs(B),18,'left-msb')';
                %此时A_bi为0,1数组
                A_bi(1) = (A<0);%A<0意味着第一位为负
                A_gray = bi2graycode(A_bi);
                B_gray = bi2graycode(B_bi);
                gene(i,:,:) = [A_gray,B_gray];
            end
            self.gene = gene;
        end
        
        function param = gene2param(self)
            len = size(self.gene,1);
            param = zeros(len,2);
            for i = 1:len
                A = biarr2de(gray2biarr(self.gene(i,:,1)));
                B = biarr2de(gray2biarr(self.gene(i,:,2)));
                param(i,:) = [A,B]*1e-6;
            end
            self.param = param;
        end
        
        %计算params的fitnessfun
        function  fitness = fitnessfun(self,S,fitfun)
            params = self.param;
            wl = S.wl;
            N = S.N;
            t = S.t;
            wh = zeros(1,length(params));
            th = zeros(1,length(params));
            for i = 1:length(params)
                %B<0适应度为0
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
            fitness = fitfun.Get_fit(S1,S);
        end
        
    end
end

%把二进制列数组转化为整数,第一位表示正负
function de = biarr2de(biarr)
    
    de = bi2de(biarr(2:end)','left-msb');
    if biarr(1)==1
        de = -de;
    end
end

%输入bi_str为列数组,返回gray码数组
function gray_code = bi2graycode(bi_arr)
    gray_code = false(length(bi_arr),1);
%     bi_code = arrayfun(@(x) str2double(x),bi_str);
    gray_code(1) = bi_arr(1);
    %第一位保留,其他的对前一位和这一位作异或运算
    for i = 2:length(bi_arr)
        gray_code(i) = xor(bi_arr(i-1),bi_arr(i));
    end
end

%输入为gray码列数组，返回2进制数组
function bi_arr = gray2biarr(gray_arr)
    len = length(gray_arr);
    bi_arr = false(len,1);
    bi_arr(1) = gray_arr(1);
    for i = 2:len
        bi_arr(i) = xor(bi_arr(i-1),gray_arr(i));
    end
end
