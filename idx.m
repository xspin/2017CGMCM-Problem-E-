% J1-62, D1-2, F01-70
% 1-62, 63-64, 65-134
function i = idx(s)
    if isnumeric(s)
        if s < 63
            i = sprintf('J%02d', s);
        elseif s < 65
            i = sprintf('D%02d', s-62);
        elseif s < 71
            i = sprintf('Z%02d', s-64);
        else
            i = sprintf('F%02d', s-70);
        end
        return
    end
    if s(1) == 'J'
        i = 0;
    elseif s(1)=='D'
        i = 62;
    elseif s(1)=='Z'
        i = 64;
    elseif s(1)=='F'
        i = 70;
    else
        return
    end
    i = i + str2num(s(2:end));
end