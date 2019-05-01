%�Ŵ��㷨
%����ȷ��4λС���Ĳ���A,B��10000ȡ����ת��Ϊ�����ƴ�
%��࿼��6λ����,ʹ��gray���봢��
%ÿһ���Ļ���洢��һ��3άlogical����N_c*18*2��
%����Ҳ�洢��2ά�����У�N_c*2
%��ʼ�趨N����(ǰ������ͨ��Ѱ���㷨ֻ��һ������,���������������)
%gray��ʹ�����鷽��תΪ�����ƴ�
%����Ⱦɫ���Ϊlogical array
classdef Individual<matlab.mixin.Copyable
    properties
        gene   %����(�ֿ�������)
        param        %����,����,N_c*2
        chromosome  %Ⱦɫ��,Ϊ��gene���ӳ�һ���ַ���,Ϊ������
        N_c
        fitness
    end
    methods
        
        %��ʼ��
        %S��Ӧobject function
        function self = Individual(N_c,wh_max,S,fitfun)
            self.N_c = N_c;
            self.param = zeros(N_c,2);
%             wh_pre = wh_max*rand(1,N_c);
            %�˵ĳ�ʼֵ���ϼ�С
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
        
        %��Ⱦɫ����º�Ի���Ͳ������и���
        %�ڲ������������������¸���ʱ�õ�
         function Renew(self,S,fitfun)
            self.segment_chro();
            self.gene2param();
            self.fitness = self.fitnessfun(S,fitfun);
         end
         
         %�����º�Ⱦɫ��ָ�Ϊ����
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
             %ÿ���˶�Ӧ�Ļ��򳤶ȣ�ӦΪ36
             len = length(chro)/self.N_c;
             for i = 1:self.N_c
                 chro((i-1)*len+1:i*len) = [self.gene(i,:,1)',...
                     self.gene(i,:,2)'];
             end
             self.chromosome = chro;
         end
         
        %����paramΪһ����ά����(��)
        %geneΪһ��N_c*n*2�ľ��󣬴���A,B��������
        function gene = param2gene(self)
            gene = false(self.N_c,18,2);
            for i = 1:length(self.param)
                A = floor(self.param(i,1)*1e6);
                B = floor(self.param(i,2)*1e6);
                A_bi = de2bi(abs(A),18,'left-msb')';
                B_bi = de2bi(abs(B),18,'left-msb')';
                %��ʱA_biΪ0,1����
                A_bi(1) = (A<0);%A<0��ζ�ŵ�һλΪ��
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
        
        %����params��fitnessfun
        function  fitness = fitnessfun(self,S,fitfun)
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
            fitness = fitfun.Get_fit(S1,S);
        end
        
    end
end

%�Ѷ�����������ת��Ϊ����,��һλ��ʾ����
function de = biarr2de(biarr)
    
    de = bi2de(biarr(2:end)','left-msb');
    if biarr(1)==1
        de = -de;
    end
end

%����bi_strΪ������,����gray������
function gray_code = bi2graycode(bi_arr)
    gray_code = false(length(bi_arr),1);
%     bi_code = arrayfun(@(x) str2double(x),bi_str);
    gray_code(1) = bi_arr(1);
    %��һλ����,�����Ķ�ǰһλ����һλ���������
    for i = 2:length(bi_arr)
        gray_code(i) = xor(bi_arr(i-1),bi_arr(i));
    end
end

%����Ϊgray�������飬����2��������
function bi_arr = gray2biarr(gray_arr)
    len = length(gray_arr);
    bi_arr = false(len,1);
    bi_arr(1) = gray_arr(1);
    for i = 2:len
        bi_arr(i) = xor(bi_arr(i-1),gray_arr(i));
    end
end
