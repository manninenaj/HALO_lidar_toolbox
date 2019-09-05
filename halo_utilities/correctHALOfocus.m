function data = correctHALOfocus(site,DATE,abc,data)

if any(strcmp(site,{'hyytiala','sodankyla','uto','southbank'}))
    C = getconfig(site,DATE);
    [data.beta_raw,~] = correct_focus(C.(['focus_' abc]), data);
    data.beta = data.beta_raw;
    data.v = data.v_raw;
end
end

