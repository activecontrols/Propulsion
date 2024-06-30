function apx(value)
    %This is a quick function to aproximate the symbolic values; it is intended to be used for diagnosal purposes only. If
    %keepUnits is a boolean for if you want the units
    if true
        [val, unit] = separateUnits(value);
        val = double(val);
        fprintf("%d ", val);
        display(unit);
        fprintf("\n")
    end
end