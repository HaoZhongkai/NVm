function chro_out = inverse(chro)
    %��λ
    len = length(chro);
    if rand>0
        inverse_fragment = sort(randi(len,1,2));
        chro_out = chro;
        chro_out(inverse_fragment) = fliplr(chro(inverse_fragment));
    else
        %������Bͻ��
        len_gene = 36;
        down_adjust = randi(length(chro)/len_gene)-1;
        chro_out = chro;
        chro_out(down_adjust*len_gene+19:(down_adjust+1)*len_gene) = ...
            mutate(chro(down_adjust*len_gene+19:(down_adjust+1)*len_gene),...
            0.4);
    end
end