function patch = patch_data_function

patch = struct(...
    'nv', 0, ...        % patch normal vector
    'px', 0, ...        % patch ring points
    'num_grid', 0, ...  % number of grid points
    'tv', 0, ...        % patch tangential vectors
    'dtv', 0, ...       % patch grid size
    'p', 0, ...         % grid point positions
    'v', 0 ...          % velocites at grid
    ); 

end