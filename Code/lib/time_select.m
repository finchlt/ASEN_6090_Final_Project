function tk = time_select(t,toe)
  
    tk_tmp = t - toe;
    
    if tk_tmp > 302400
        tk = tk_tmp - 604800;
    elseif tk_tmp < -302400
        tk = tk_tmp + 604800;
    else
        tk = tk_tmp;
    end
end
