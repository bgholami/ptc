function [p, f] = parse_inlet_plt_file(input_file)

fid = fopen(input_file, 'r');
N = fscanf(fid, '%d', 3);
Nel = N(1); % number of elements
Npo = N(2); % number of points
Nfa = N(3); % number of facets

% skip over Nel+1 lines (first line + elements)
for i = 1:Nel+1
    fgetl(fid);
end

% read points into p(1:num_dim, Npo)
for i = 1:Npo
    index = fscanf(fid, '%d', 1);
    data  = fscanf(fid, '%f', 3);
    p(1:3, index) = data(1:3);
end

% read facets into f(1:3, Nfa), i.e. 3 points per facet
for i = 1:Nfa
    data = fscanf(fid, '%d', 5);
    f(1:3, i) = data(1:3);
end

fclose(fid);