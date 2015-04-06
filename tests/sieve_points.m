function [flist, bc_points] = sieve_points(f)

% This should return facets that form the surface. Ideally, facets could be
% indetified using IDs in the input file:
% wall = 1, inlet = 2, outlets > 2

flist = zeros(3, 0);
for i = 1:3
    bc_points{i} = []; % 1 inlet, 2:3 outlet
end
counter = 0;
for i = 1:size(f, 2)
    % if all points belong to surface
    if (facet_status(f(:, i)) == 0)
        counter = counter + 1;
        flist(:, counter) = f(:, i);
    else % i.e. inlet/outlet facet
        for ip = 1:3
            if (p_status(f(ip, i)) == 0) 
                % i.e. if this non-wall facet has a point on wall
                bc_points{facet_status(f(:, i))} = [bc_points{facet_status(f(:, i))}, f(ip, i)];
            end
        end
    end
end

% make bc_points enteies unique
for i = 1:3
    bc_points{i} = unique(bc_points{i});
end
end


function status = facet_status(plist)

for i = 1:3
    pstatus(i) = p_status(plist(i));
end

if (sum(pstatus) == 0)
    status = 0;
    return
else
    temp = pstatus(pstatus~=0);
    status = temp(1);
    return
end

end


function status = p_status(id)

% since the current input file does not have correct IDs. We rely on
% information that is not automatically generated.
% wall = 0, inlet = 1, outlets > 1
%inlet_ring = 5222:5266;
%outlet1_ring = 5378:5398;
%outlet2_ring  = 5424:5446;
inlet_ring = 5222:5377;
outlet1_ring = 5378:5423;
outlet2_ring  = 5424:5473;

status = 0;
if (id < 5222)
    status = 0;
    return
else
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