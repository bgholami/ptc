% This tests the algorithm of Sazonov2011modelling
clear
clc

% define cartiesian grid
mu = 0.04;
T  = 0.92;
a  = .23;
t = [0, 4.91254E-2, 7.81104E-2, 1.20856E-1, 2.21335E-1, 3.29924E-1, 4.21961E-1, 7.81572E-1];

% a polygon to mark the edges of the boundary
for i = 1:100
    alpha = 2*pi/100;
    bctp(i, 1) = a * cos(alpha * (i-1));
    bctp(i, 2) = a * sin(alpha * (i-1));
end


% grid
g = linspace(-a, a, 51);
dx = g(2) - g(1);
[X, Y] = meshgrid(g, g);

c = 0;
for j = 1:size(X, 1)
    for k = 1:size(X, 2)
        c = c + 1;
        vx(k, j, 1) = X(j, k);
        vx(k, j, 2) = Y(j, k);
    end
end

bc_normal_velocity = compute_Womersley_profiles(vx, bctp(:, 1:2), t, T, mu, '/home/jo/ptc/tests/fig18.csv');

% calculate the inlet velocity profile at every instant
figure
img = imread('saz.png');
imagesc([-1 1], [0 140], flipdim(img,1));
set(gca,'YDir','normal')

hold on

%color= 'cbmygrmk';
color= 'gyrrmbck';

for it = 1:length(t)
% reshape data
u = reshape(squeeze(bc_normal_velocity(it, :, :))', size(vx, 1)*size(vx, 2), 1);
for j = 1:2
    ux(:, j) = reshape(vx(:, :, j)', size(vx, 1)*size(vx, 2), 1);
end

sgn = sign(ux(:, 1));
sgn(sgn==0) = 1;
cr = sgn.*(ux(:, 1).^2 + ux(:, 2).^2).^0.5;
cr = cr / a;
plot(cr, u, 'r.')
end
