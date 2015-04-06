function status = pp_status(id)

% since the current input file does not have correct IDs. We rely on
% information that is not automatically generated.
% wall = 0, inlet = 1, outlets > 1
inlet_ring = 5222:5266;
outlet1_ring = 5378:5398;
outlet2_ring  = 5424:5446;

status = 0;
if (id < 5222)
    status = 0;
    return
else
    id
    if (min(abs(inlet_ring - id)) == 0)
        status = 1;
        return
    elseif (min(abs(outlet1_ring - id)) == 0)
        status = 2;
        return
    elseif (min(abs(outlet2_ring - id)) == 0)
        status = 3;
        return
    end    
end

end