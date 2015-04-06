function test_surface(f, p)

pp = 0;
counter = 0;
for i = 1:size(p, 2)
    if min(min(abs(f - i))) == 0
        counter = counter + 1;
        pp(counter) = i;
    end
end

fid = fopen('test.csv', 'w');
for i = 1:length(pp)
    fprintf(fid, '%f %f %f\n', p(1:3, pp(i)));
end
fclose(fid);

inlet_ring = 5222:5266;
outlet1_ring = 5378:5398;
outlet2_ring  = 5424:5446;

scatter3(p(1, pp(inlet_ring)), p(2, pp(inlet_ring)), p(3, pp(inlet_ring)), '.');
hold on;
scatter3(p(1, pp(outlet1_ring)), p(2, pp(outlet1_ring)), p(3, pp(outlet1_ring)), '.');
scatter3(p(1, pp(outlet2_ring)), p(2, pp(outlet2_ring)), p(3, pp(outlet2_ring)), '.');
axis equal