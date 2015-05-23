function bc_velocity = compute_Womersley_profiles(pp, bcp, time, T, mu, input_file)
% This follows the algorithm of Sazonov2011modelling to calculate the
% transient velocity profiles at boundaries

% read the velocity waveform
velocity_waveform = load(input_file); % max velocity over time

% define cartiesian grid
N = 24; % number of modes? TODO

% compute complex amplitudes
U = compute_complex_amplitudes(velocity_waveform, 1000);

% reshape data
for j = 1:2
    tvx(:, j) = reshape(pp(:, :, j), size(pp, 1)*size(pp, 2), 1);
end

% solve the boundary value problem
[v, vx, mif] = solve_boundary_value_problem(tvx, bcp, N, T, mu, size(pp, 2));

% calculate the inlet velocity profile at every instant
for it = 1:length(time)
    u = ones(length(vx),1) * U(1) .* v(:, 1);
    for n = 1:N
        wn = 2*pi*n/T;
        u = u + 2 * real(U(n+1) * v(:, n+1) * exp(1i * wn * time(it)));
    end
    %sum(isnan(u))
    tbcv = zeros(length(tvx), 1);
    tbcv(mif) = u;
    bc_velocity(it, :, :) = reshape(tbcv, size(pp, 1), size(pp, 2));
    %sgn = sign(vx(:, 1));
    %sgn(sgn==0) = 1;
    %cr = sgn.*(vx(:, 1).^2 + vx(:, 2).^2).^0.5;
    %cr = cr / a;
    %ppp = polyfit(cr, u, 9);
    %cr = sort(cr);
    %pvv = polyval(ppp, cr);
    %plot(cr, pvv, color(it),'linewidth', 2);
    %fl(it) = compute_flow_rate3D(cr, pvv);
    %plot(cr, pvv, '--k','linewidth', 2);
end
