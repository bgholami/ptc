function write_bcdef_file(patch, num_inflow, num_outflow, id, time, outfile)

num_inout = num_inflow + num_outflow;
num_time = length(time);

fid = fopen(outfile, 'w');
fprintf(fid,'%d %d\n', num_inflow, num_outflow);
fprintf(fid,'%d %d %d\n', size(patch(1).px, 1), size(patch(2).px, 1), size(patch(3).px, 1)); % write all patches in one line
fprintf(fid,'%d ', id);
fprintf(fid,'\n');

for i = 1:num_inout
    fprintf(fid,'%15.16f %15.16f %15.16f\n', patch(i).nv);
    for ip = 1:size(patch(i).px, 1)
        fprintf(fid,'%15.16f %15.16f %15.16f\n', patch(i).px(ip, :));
    end
end

max_num_grid1 = -1;
max_num_grid2 = -1;
for i = 1:num_inout
    max_num_grid1 = max(max_num_grid1, patch(i).num_grid1);
    max_num_grid2 = max(max_num_grid2, patch(i).num_grid2);
end

fprintf(fid,'%d %d %d\n', num_time, max_num_grid1, max_num_grid2);
fprintf(fid,'%15.16f\n', time);

for i = 1:num_inout
    
    % write
    fprintf(fid,'%d %d\n', patch(i).num_grid1, patch(i).num_grid2);
    fprintf(fid,'%15.16f %15.16f %15.16f\n', patch(i).tv1);
    fprintf(fid,'%15.16f %15.16f %15.16f\n', patch(i).tv2);
    fprintf(fid,'%15.16f %15.16f\n', patch(i).dtv1, patch(i).dtv2);
    fprintf(fid,'%15.16f %15.16f %15.16f\n', patch(i).p(1,1,:));
    
    for it = 1:num_time
        for ii = 1:patch(i).num_grid1
            for jj = 1:patch(i).num_grid2
                fprintf(fid,'%15.16f %15.16f %15.16f\n', patch(i).v(1:3, it, ii, jj));
            end
        end
    end
    
end

fclose(fid);

end