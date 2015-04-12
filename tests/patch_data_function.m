function patch = patch_data_function

patch = struct(...
    'nv', 0, ...        % patch normal vector    
    'tv1', 0, ...       % first patch tangential vector
    'tv2', 0, ...       % second patch tangential vector
    'px', 0, ...        % patch ring points
    'num_grid1', 0, ... % number of grid points (along tv1)
    'num_grid2', 0, ... % number of grid points (along tv2)
    'dtv1', 0, ...      % patch grid size (along tv1)
    'dtv2', 0, ...      % patch grid size (along tv2)
    'p', 0, ...         % grid point positions
    'v', 0 ...          % velocites at grid
    ); 

end