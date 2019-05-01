function de = bistr2de(bistr)
    ref = arrayfun(@(x) str2double(x),bistr);
    de = bi2de(ref(2:end)','left-msb');
    if bistr(1)==1
        de = -de;
    end
end

