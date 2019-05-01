%输入为1*n logical array，输出相同
function chro_out = mutate(chro,mutate_ratio)
    len = length(chro);
    Mutate0 = rand(1,len);
    chro_out = chro;
    for i = 1:len
        if Mutate0(i)<mutate_ratio
            chro_out(i) = 1-chro(i);
        end
    end
end