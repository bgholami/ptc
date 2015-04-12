% This follows the algorithm of Sazonov2011modelling to calculate the
% transient velocity profiles at boundaries
clear
clc

% read the velocity waveform
velocity_waveform = load('/home/jo/ptc/tests/womer/fig18.csv'); % max velocity over time

% define cartiesian grid
mu = 0.04;
T  = 0.92;
a  = .23;
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
figure
img = imread('saz.png');
imagesc([-1 1], [0 140], flipdim(img,1));
set(gca,'YDir','normal')

hold on
t = [0, 4.91254E-2, 7.81104E-2, 1.20856E-1, 2.21335E-1, 3.29924E-1, 4.21961E-1, 7.81572E-1];
%color= 'cbmygrmk';
color= 'gyrrmbck';

for it = 1:length(t)
u = ones(length(vx),1) * U(1) .* v(:, 1);
for n = 1:N
    wn = 2*pi*n/T;
    u = u + 2 * real(U(n+1) * v(:, n+1) * exp(1i * wn * t(it)));
end

%figure
%scatter3(vx(:, 1), vx(:, 2), u, 'r.')
sgn = sign(vx(:, 1));
sgn(sgn==0) = 1;
cr = sgn.*(vx(:, 1).^2 + vx(:, 2).^2).^0.5;
cr = cr / 0.23;
if (0)
    plot(cr, u, 'r.')
else
    ppp = polyfit(cr, u, 9);
    cr = sort(cr);
    pvv = polyval(ppp, cr);
    %plot(cr, pvv, color(it),'linewidth', 2);
    plot(cr, pvv, '--k','linewidth', 2);
end
end
