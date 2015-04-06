clear
clc

% open input plt file
[p, f] = parse_inlet_plt_file('Carotid01Rm_427_3D.plt');
Npo = size(p, 2); % number of points
Nfa = size(f, 2); % number of facets

% optional: test surface
% test_surface(f, p);

% indetify surface and inlet/outlet points
[flist, bc_points] = sieve_points(f);

%%% create geometry input     
% sort according to this coordinate
sort_cor = 3;
near_wall = 6.0844100E-001; % mm
create_geometry_input(p, flist, sort_cor, near_wall);

% create inlet/outlet input
%create_bcdef_input(bc_points)
    
for i = 1:3
    scatter3(p(1, bc_points{i}), p(2, bc_points{i}), p(3, bc_points{i}), '.')
    hold on
end