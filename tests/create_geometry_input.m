function create_geometry_input(p, flist, sort_cor, dnear_wall, outfile)

Nfa = size(flist, 2);
Npo = length(unique(flist));
nv = 0;

facet_vlist = zeros(Nfa,3);
vertex_flist = zeros(Npo,20);
v = zeros(Npo,3);

for fi = 1:size(flist, 2)
    
    % calculate facet normal
    fnormal(fi, :) = calculate_facet_normal(p(:, flist(1:3, fi)));
    
    for ip = 1:3        
        %vt = str2num(tl(8:end));
        vt = p(:, flist(ip, fi));
        
        % create list data
        for i = 1:3
            temp(:, i) = v(:, i) - vt(i);
        end
        [minv, minl] = min(sum(abs(temp), 2));
        if (minv == 0)
            % means this vertex has already been stored
            facet_vlist(fi, :) = [minl, facet_vlist(fi, 1:2)];
            vertex_flist(minl, :) = [fi, vertex_flist(minl, 1:end-1)];
        else
            % means this is a new vertex
            nv = nv + 1;
            v(nv, :) = vt;
            facet_vlist(fi, :) = [nv, facet_vlist(fi, 1:2)];
            vertex_flist(nv, :) = [fi, vertex_flist(nv, 1:end-1)];
        end        
    end
end

% sort
if (1)
      
    [~, ix] = sort(v(:, sort_cor));
    for i = 1:length(ix)
        rix(ix(i)) = i;
    end
    
    sv = v(ix, :);
    svertex_flist = vertex_flist(ix, :);
    for i = 1:size(facet_vlist, 1)
        for j = 1:size(facet_vlist, 2)
            sfacet_vlist(i, j) = rix(facet_vlist(i, j));
        end
    end
    sfnormal = fnormal;
    
    % find the smallest dx larger than dnear_wall that suits the geometry
    %dnear_wall = 6.0844100E-001; % mm
    dx = dnear_wall * 1.1;
    
    xmin = min(sv(:, sort_cor));
    xmax = max(sv(:, sort_cor));
    
    % using dx calculate points of division
    divs = zeros(floor((xmax - xmin) / dx) + 2, 1);
    for i = 1:length(sv)
        ind = floor((sv(i, sort_cor) - xmin) / dx) + 2;
        divs(ind) = divs(ind) + 1; % it will work because sv is sorted
    end
    for i = 1:length(divs)-1
        divs(i+1) = divs(i) + divs(i+1);
    end
    
    
else
    
    sv = v;
    svertex_flist = vertex_flist;
    sfacet_vlist = facet_vlist;
    sfnormal = fnormal;
    
end

% test the grid
if (1)
    
    v_inc = zeros(length(sv), 1);
    
    for i = 1:size(sfacet_vlist, 1)
        res = 0;
        for j = 1:size(sfacet_vlist, 2)
            vi = sfacet_vlist(i, j);
            v_inc(vi) = 1; % every vertex should belong to at least one facet
            res = res + min(abs((svertex_flist(vi, :) - i)));
        end
        if (res ~= 0)
            error('inconsistent grid at facets');
        end
    end
    if (sum(v_inc) ~= length(sv))
        error('inconsistent grid at vertices');
    end
    
end

%write output
if (1)
    outfile1 = outfile;
    outfile2 = [outfile, '.csv'];
    write_arbitrary_data_sorted(outfile1, outfile2, sv, svertex_flist, sfacet_vlist, sfnormal, dx, divs, sort_cor)
end


end


function fnormal = calculate_facet_normal(points)

v1 = points(:, 1) - points(:, 2);
v2 = points(:, 1) - points(:, 3);
vn = cross(v1, v2);
fnormal = vn / norm(vn);

end