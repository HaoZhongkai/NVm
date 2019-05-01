function [chro_out1,chro_out2] = Segment_crossover(chro1,chro2)
    %单基因长度36
    %交换一个单基因B
    len = 36;
    segment1 = randi(length(chro1)/len-1);
    segment1 = segment1*len+1:(segment1+1)*len;
    segment2 = randi(length(chro1)/len-1);
    segment2 = segment2*len+1:(segment2+1)*len;
    chro_out1 = chro1;
    chro_out2 = chro2;
    chro_out1(segment1) = chro2(segment2);
    chro_out2(segment2) = chro1(segment1);
end