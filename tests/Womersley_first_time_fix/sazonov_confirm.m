clear
clc
%load matlab.mat
mu = 0.04;
T  = 0.92;
a = .23;

N = 24;

r = .23;
c = 1;
x(1,1) = 0;
y(1,1) = 0;
dx = 0.1 * .23;

na1 = floor(r/dx) + 1;
da1 = r / na1;

for j = 1:na1
    ra = da1 * j;
    na2 = floor(2*pi*ra/dx) + 1;
    da2 = 2*pi / na2;
    
    for k = 1:na2
        c = c + 1;
        
        a = da2 * (k-1);
        x(c,1) = ra * cos(a);
        y(c,1) = ra * sin(a);
        
    end
end

DT = delaunayTriangulation([x,y]);
fe = freeBoundary(DT);
fe = fe(:, 1);
vx = DT.Points;
fvlist = DT.ConnectivityList;

triplot(DT);
hold on;
plot(x(fe),y(fe),'-r','LineWidth',2) ; 
axis equal
hold off;
a = .23;
for n = 0:N
    
    if (n == 0)
         
    for k = 1:length(vx)
            
            rr = (vx(k, 1)^2 + vx(k, 2)^2)^0.5;
            v(k, n+1) = 0;
            v(k, n+1) = (1 - rr^2 / a^2);
            
    end
        
    else
    
    wn = 2*pi*n/T;
    kn = (wn/2/mu)^0.5 * (1-1i);
    %kn = (-1i * wn / mu)^0.5;
    
    for k = 1:length(vx)
            
            rr = (vx(k, 1)^2 + vx(k, 2)^2)^0.5;
            v(k, n+1) = 0;
            v(k, n+1) = (1 - besselj(0, kn * rr) / besselj(0, kn * a)) / (1 - 1/besselj(0, kn * a));
            
    end
    
    end
end

% complex amplitudes
if (0) % data provided in Sazonov2011
    someu = load('/home/gholami/work/tests/sazonov_tableCI');
    U = someu(:, 3) .* exp(1i*someu(:, 4));
else % calculate the complex amplitude data
    b = load('/home/gholami/work/tests/fig18.csv'); % max velocity over time
    U = compute_complex_amplitudes(b, 100);
end

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
%plot(cr, u, 'r.')
ppp = polyfit(cr, u, 9);
cr = sort(cr);
pvv = polyval(ppp, cr);
%plot(cr, pvv, color(it),'linewidth', 2);
plot(cr, pvv, '--k','linewidth', 2);
end