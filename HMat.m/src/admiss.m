function ad = admiss(x,y,type_admiss)

if (type_admiss == 'S')
    if max(x-y) > 1
        ad = true;
    else
        ad = false;
    end
elseif type_admiss == 'E'
    if sum(x~=y) > 1
        ad = true;
    else
        ad = false;
    end
elseif type_admiss == 'W'
    if sum(x~=y) >= 1
        ad = true;
    else
        ad = false;
    end
end

end