%输入都是行向量字符串
function [chro_out1,chro_out2] = uniform_crossover(chro1,chro2)
    len = length(chro1);
    Crossover0 = rand(1,len);
    chro_out1 = chro1;
    chro_out2 = chro2;
    for i = 1:len
        if Crossover0(i)<1/2
            chro_out1(i) = chro2(i);
            chro_out2(i) = chro1(i);
        end
    end
end