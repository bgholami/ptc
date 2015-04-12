% This follows the algorithm of Sazonov2011modelling to calculate the
% transient velocity profiles at boundaries
clear
clc

% read the velocity waveform
velocity_waveform = load('hinds_vmax.data'); % max velocity over time

Wo = 3.0;
mu = 0.008;
a  = 0.635/2;
w = mu * (Wo / a)^2;
T = 1 / w;
%T = velocity_waveform(end, 1) - velocity_waveform(1, 1);
velocity_waveform(:, 1) = (velocity_waveform(:, 1) - velocity_waveform(1, 1)) / ...
    (velocity_waveform(end, 1) - velocity_waveform(1, 1)) * T;
Wo = ((1/T) * a^2 / mu)^0.5

flowrate = velocity_waveform(:, 2)*2/3*pi*a*a;
% define cartiesian grid
N = 24;
nn = 51;
g = linspace(-a, a, nn);
dx = g(2) - g(1);
[X, Y] = meshgrid(g, g);

c = 0;
for j = 1:nn
    for k = 1:nn
        c = c + 1;
        vx(c, 1) = X(j, k);
        vx(c, 2) = Y(j, k);
    end
end

% compute complex amplitudes
U = compute_complex_amplitudes(velocity_waveform, 1000);

% solve the boundary value problem
[v, vx] = solve_boundary_value_problem(vx, a, N, T, mu, nn, dx);

% calculate the inlet velocity profile at every instant
%t = [0, 1.78, 3.724, 5.368, 6.879, 8.104];
t = [0, 0.2725, 0.5701, 0.8218, 1.0531, 1.2406];
%t = linspace(velocity_waveform(1, 1), velocity_waveform(end, 1), 7);
color= 'rbgkmykkk';
hold on
for it = 1:length(t)
    u = ones(length(vx),1) * U(1) .* v(:, 1);
    for n = 1:N
        wn = 2*pi*n/T;
        u = u + 2 * real(U(n+1) * v(:, n+1) * exp(1i * wn * t(it)));
    end
    
    sgn = sign(vx(:, 1));
    sgn(sgn==0) = 1;
    cr = sgn.*(vx(:, 1).^2 + vx(:, 2).^2).^0.5;
    %cr = cr / a;
    ppp = polyfit(cr, u, 9);
    cr = sort(cr);
    pvv = polyval(ppp, cr);
    plot(cr, pvv, color(it),'linewidth', 2);
    fl(it) = compute_flow_rate3D(cr, pvv);
    %plot(cr, pvv, '--k','linewidth', 2);
end

figure
plot(velocity_waveform(:, 1), flowrate)
hold on 
plot(t, fl, 'ro')
