clear
clc

WWo = [0.1, 1.0, 2.0, 5.0];
 prepare_figure(16/3, 1500) % arguments: ration, width (optinal)
for iii = 1:4
    subplot(1,4,iii, 'FontSize',15)
mu = 0.04;
T  = 0.92;
a = .23;
dt = 0.01;

x = 0:dt:T;
y = sin(pi*x) + 1;
velocity_waveform = [x; y]';

Wo = WWo(iii);
w = 1 / T;
mu = w / (Wo / a)^2;
Wo = (w * a^2 / mu)^0.5

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
U = compute_complex_amplitudes(velocity_waveform, 100);

% solve the boundary value problem
[v, vx] = solve_boundary_value_problem(vx, a, N, T, mu, nn, dx);

% calculate the inlet velocity profile at every instant
t = linspace(0, T, 6);
color= 'krbgkmykkk';
hold on
for it = 2:length(t)
    u = ones(length(vx),1) * U(1) .* v(:, 1);
    for n = 1:N
        wn = 2*pi*n/T;
        u = u + 2 * real(U(n+1) * v(:, n+1) * exp(1i * wn * t(it)));
    end
    
    sgn = sign(vx(:, 1));
    sgn(sgn==0) = 1;
    cr = sgn.*(vx(:, 1).^2 + vx(:, 2).^2).^0.5;
    cr = cr / a;
    ppp = polyfit(cr, u, 9);
    cr = sort(cr);
    pvv = polyval(ppp, cr);
    plot(cr, pvv, color(it),'linewidth', 2);
    %plot(cr, pvv, '--k','linewidth', 2);
end
box on 
grid on 
xlabel('r')
ylabel('nomalized velocity')
title(['Wo = ', num2str(Wo)])
ylim([-0.5 2.5])
end

