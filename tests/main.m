clear
clc

%%%%% input specification    

% sort according to this coordinate
sort_cor = 3;
near_wall = 6.0844100E-001;

%%%  files
input_CCA_geometry_plt = 'Carotid01Rm_427_3D.plt';
% velocity waveforms
input_inlet_vw = 'CCA_pulse/cca.csv';
input_outlet1_vw = 'CCA_pulse/ica.csv';
input_outlet2_vw = 'CCA_pulse/eca.csv';

output_bc_definition = 'bcdef_CCA01.data';
output_geo_definition = 'geom_CCA01_sorted.data';

% characteristic dx to determine number of patch grids
T = 0.8896; % waveform period
mu = 0.008; % viscosity

% time-varying velocity data
min_num_grid = 11; % at least this many grid points in each direction
%time = linspace(0, 9, 5); % time instances to write velocity at patches 
time = [0,0.0553237500000000,0.110647500000000,0.165971250000000,0.221295000000000,0.276618750000000,0.331942500000000,0.387266250000000,0.442590000000000,0.497913750000000,0.553237500000000,0.608561250000000,0.663885000000000,0.719208750000000,0.774532500000000,0.829856250000000,0.885180000000000];

%%%------------------------------------------


if (1)
    % open input plt file
    [p, f] = parse_inlet_plt_file(input_CCA_geometry_plt);
    
    % optional: test surface
    % test_surface(f, p);
    
    % indetify surface and inlet/outlet points
    [flist, bc_points] = sieve_points(f);
    
    % create geometry input
    create_geometry_input(p, flist, sort_cor, near_wall, output_geo_definition);
else
    load workspace.mat
end

% create inlet/outlet input
create_bcdef_input(p, bc_points, time, T, mu, min_num_grid, output_bc_definition, {input_inlet_vw, input_outlet1_vw, input_outlet2_vw});

for i = 1:3
    scatter3(p(1, bc_points{i}), p(2, bc_points{i}), p(3, bc_points{i}), '.')
    hold on
end