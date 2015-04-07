function create_bcdef_input(p, bc_points, outfile)

num_inflow = 1;
num_outflow = 2;

% characteristic dx to determine number of patch grids
char_dx = 0.09170 / 3.0;

% time-varying velocity data
num_grid_1(1:3) = [15 11 9];
num_grid_2(1:3) = [13 11 9];
num_time   = 5;
time = linspace(0, 9, num_time);

% TODO: 3 input velocities, one per patch, each in this form
%v(1:3, 1:num_time, 1:num_grid_1(i), 1:num_grid_2(i))

% ---------------------------------------------------

num_inout = num_inflow + num_outflow;

% one patch for each inflow/outflow boundary
patch(1:num_inout) = patch_data_function;

% ID of each patch (abs(id) = 1:num_inflow+num_outflow, sign(id) = 1 for inlet, -1 for outlet)
id = 1:num_inflow+num_outflow;
id(num_inflow+1:end) = -1 * id(num_inflow+1:end);

% midpoint in domain
mp = sum(p, 2) / length(p);

fid = fopen(outfile, 'w');
fprintf(fid,'%d %d\n', num_inflow, num_outflow);
fprintf(fid,'%d %d %d\n', length(bc_points{1}), length(bc_points{2}), length(bc_points{3})); % write all patches in one line
fprintf(fid,'%d ', id);
fprintf(fid,'\n');

% boundary ring points
for i = 1:num_inout
    patch(i).px = p(1:3, bc_points{i})';
    
    % boundary normals
    dn = calculate_bc_patch_normal(patch(i).px, mp);
    
    patch(i).nv = dn;
    
    
    fprintf(fid,'%15.16f %15.16f %15.16f\n', patch(i).nv);
    for ip = 1:size(patch(i).px, 1)
        fprintf(fid,'%15.16f %15.16f %15.16f\n', patch(i).px(ip, :));
    end
    
end

fprintf(fid,'%d %d %d\n', num_time, max(num_grid_1), max(num_grid_2));
fprintf(fid,'%15.16f\n', time);

for i = 1:num_inout
    
    nv = patch(i).nv;
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
    
    % define the patch in this coordinate system: t1, t2, nv
    for ii = 1:length(patch(i).px)
        qx1(ii) = dot(patch(i).px(ii, :), t1);
        qx2(ii) = dot(patch(i).px(ii, :), t2);
    end
    lim1 = [min(qx1), max(qx1)];
    lim2 = [min(qx2), max(qx2)];
    
    x1 = linspace(lim1(1), lim1(2), num_grid_1(i));
    x2 = linspace(lim2(1), lim2(2), num_grid_2(i));
    
    dt1 = x1(2) - x1(1);
    dt2 = x2(2) - x2(1);
    
    
    for ii = 1:num_grid_1(i)
        for jj = 1:num_grid_2(i)
            pp(ii, jj, 1:3) = x1(ii) * t1 + x2(jj) * t2;
        end
    end
    
    % TODO: replace this line with input
    v(1:3, 1:num_time, 1:num_grid_1(i), 1:num_grid_2(i)) = 0.0;
    
    % write
    fprintf(fid,'%d %d\n', num_grid_1(i), num_grid_2(i));
    fprintf(fid,'%15.16f %15.16f %15.16f\n', t1);
    fprintf(fid,'%15.16f %15.16f %15.16f\n', t2);
    fprintf(fid,'%15.16f %15.16f\n', dt1, dt2);
    fprintf(fid,'%15.16f %15.16f %15.16f\n', pp(1,1,:));
    
    for it = 1:num_time
        for ii = 1:num_grid_1(i)
            for jj = 1:num_grid_2(i)
                fprintf(fid,'%15.16f %15.16f %15.16f\n', v(1:3, it, ii, jj));
            end
        end
    end
    
end


fclose(fid);

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