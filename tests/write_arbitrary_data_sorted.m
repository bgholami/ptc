function write_arbitrary_data_sorted(outfile1, outfile2, v, vertex_flist, facet_vlist, fnormal, dx, divs, sort_cor)

% maximum degree of a vertex (max number of facets for a vertex)
max_deg = sum(min(abs(max(vertex_flist,[],1)), 1));

fid1 = fopen(outfile1, 'w');
fid2 = fopen(outfile2, 'w');

fprintf(fid1,'%d %d\n', size(v, 1), max_deg);

for i = 1:size(v, 1)
    fprintf(fid1,'%15.16f %15.16f %15.16f\n', v(i, :));
    fprintf(fid2,'%15.16f %15.16f %15.16f\n', v(i, :));
end

format = ['%d'];
for i = 1:max_deg-1
    format = [format ' %d'];
end
format = [format '\n'];
for i = 1:size(v, 1)
    fprintf(fid1, format, vertex_flist(i, 1:max_deg));
end

fprintf(fid1,'%d\n', size(facet_vlist, 1));
for i = 1:size(facet_vlist, 1)
    fprintf(fid1,'%d %d %d\n', facet_vlist(i, 1:3));
end

for i = 1:size(facet_vlist, 1)
    fprintf(fid1,'%15.16f %15.16f %15.16f\n', fnormal(i, :));
end

fprintf(fid1,'%15.16f\n', dx);
fprintf(fid1,'%d %d\n', size(divs, 1), sort_cor);
for i = 1:size(divs, 1)
    fprintf(fid1,'%d\n', divs(i));
end


fclose(fid1);
fclose(fid2);

end