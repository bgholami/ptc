function create_bcdef_input(p, bc_points, time, T, mu, min_num_grid, outfile, input_file_list)

num_inflow = 1;
num_outflow = 2;
num_inout = num_inflow + num_outflow;

% one patch for each inflow/outflow boundary
patch(1:num_inout) = patch_data_function;

% ID of each patch (abs(id) = 1:num_inflow+num_outflow, sign(id) = 1 for inlet, -1 for outlet)
id = 1:num_inflow+num_outflow;
id(num_inflow+1:end) = -1 * id(num_inflow+1:end);

% midpoint in domain
mp = sum(p, 2) / length(p);

% boundary ring points
for i = 1:num_inout
    patch(i).px = p(1:3, bc_points{i})';
    
    % boundary normals
    nv = calculate_bc_patch_normal(patch(i).px, mp);
    
    % one/two vectors in tangent plane (t1, t2)
    % (we assume normal vectors are not 0)
    it = 1;
    while ((it <= 3) && (nv(it) == 0.0))
        it = it + 1;
    end
    t1 = ones(1,3);
    t1(it) = (-sum(nv) + nv(it)) / nv(it);
    t1 = t1 / (dot(t1, t1) ^ 0.50);
    t2 = cross(nv, t1);
    
    patch(i).nv = nv;
    patch(i).tv1 = t1;
    patch(i).tv2 = t2;
    
    % coordinate transformation matrix
    A = [t1; t2; nv];
    
    % define the patch in this coordinate system: t1, t2, nv
    clear bctp
    for ii = 1:length(patch(i).px)
        bctp(ii, :) = patch(i).px(ii, :) / A;
    end
    
    pre_lim1 = [min(bctp(:, 1)), max(bctp(:, 1))];
    pre_lim2 = [min(bctp(:, 2)), max(bctp(:, 2))];
    
    dx = min(pre_lim1(2)-pre_lim1(1), pre_lim2(2)-pre_lim2(1)) / (min_num_grid-1);
    
    [num_grid1, lim1] = fix_grid_range(pre_lim1, dx);
    [num_grid2, lim2] = fix_grid_range(pre_lim2, dx);
    
    x1 = linspace(lim1(1), lim1(2), num_grid1);
    x2 = linspace(lim2(1), lim2(2), num_grid2);
    
    patch(i).num_grid1 = num_grid1;
    patch(i).num_grid2 = num_grid2;
    
    dt1 = x1(2) - x1(1);
    dt2 = x2(2) - x2(1);
    
    patch(i).dtv1 = dt1;
    patch(i).dtv2 = dt2;
    
    % create each point and transform back to the original coordinate system.
    for ii = 1:num_grid1
        for jj = 1:num_grid2
            % in original coordinate system
            pp(ii, jj, 1:3) = [x1(ii), x2(jj), bctp(ii, 3)] * A;
            % in transformed coordinate system
            pptp(ii, jj, 1:2) =  [x1(ii), x2(jj)];
            scatter3(pp(ii, jj, 1), pp(ii, jj, 2), pp(ii, jj, 3), 'r.')
            hold on
        end
    end
    
    patch(i).p = pp;
    % velocities at patch over time
    bc_normal_velocity = compute_Womersley_profiles(pptp(:, :, 1:2), bctp(:, 1:2), time, T, mu, input_file_list{i});

    v1 = 100;
    v2 = 100;
    if i == 1
        bc_normal_velocity(:,:,:) = v1;
    elseif i == 2
        bc_normal_velocity(:,:,:) = v2;
    else
        bc_normal_velocity(:,:,:) = (v1 * 0.5700^2 - v2 * 0.3296^2) / (0.4104^2);
    end
    
    max(max(max(bc_normal_velocity)))
    for it = 1:length(time)
        for ii = 1:num_grid1
            for jj = 1:num_grid2
                patch(i).v(1:3, it, ii, jj) = bc_normal_velocity(it, ii, jj) * nv * sign(id(i));
            end
        end
    end
    
end

write_bcdef_file(patch, num_inflow, num_outflow, id, time, outfile);


end


function nv = calculate_bc_patch_normal(px, mp)

er = 1.0;
old_dn = zeros(1, 3);

while (er > 0.0001)
    tp = randi(length(px), 1, 3);
    d1 = px(tp(1), :) - px(tp(2), :);
    d2 = px(tp(1), :) - px(tp(3), :);
    dn = cross(d2, d1);
    dn_norm = sqrt(dot(dn, dn));
    if (dn_norm > 0.0)
        dn = dn / sqrt(dot(dn, dn));
        er = norm(dn - old_dn) / dn_norm;
        old_dn = dn;
    else
        er = 1.0;
    end
end

% determine the direction (inside the domain)
% i.e. towards the midpoint
patch_mp = sum(px, 1) / length(px);
temp_v = mp' - patch_mp;
nv = dn * sign(dot(dn, temp_v));

end


function [num_grid_out, lim_out] = fix_grid_range(flim, fdx)

fnum_grid = (flim(2) - flim(1)) / fdx + 1;
% if (rem(fnum_grid, 1) ~= 0)
%     dlim = (ceil(fnum_grid) - fnum_grid) * fdx / 2.0;
%     flim(1) = flim(1) - dlim;
%     flim(2) = flim(2) + dlim;
%     fnum_grid = ceil(fnum_grid);
% end
num_grid_out = round(fnum_grid);
lim_out = flim;

end