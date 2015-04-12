clear
clc

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


%%% mapping indices between vx(full vertices) and vint(only interior)
c = 0;
mfi = -1 * ones(length(vx), 1);
for k = 1:length(vx)
    
    if ((vx(k,1)^2+vx(k,2)^2) < a.^2)
        
        c = c + 1;
        mfi(k) = c; % map full to interior
        mif(c) = k; % map interior to full
    end
    
end

v = zeros(length(vx), N+1);
for n = 0:N
    
    wn = 2*pi*n/T;
    kn2 = -1i * wn / mu;
    %kn2 = kn2 * 4.5;
    
    A = zeros(max(mfi), max(mfi));
    b = zeros(max(mfi), 1);
    
    
    for ind1 = 1:length(vx)
        
        if (mfi(ind1) > 0) % i.e. if interior
            
            % convert to JJ, KK
            JJ = floor((ind1-1)/nn) + 1;
            KK = ind1 - (JJ-1) * nn;
            
            A(mfi(ind1),mfi(ind1)) = -4 + dx*dx*kn2;
            b(mfi(ind1)) = -dx*dx;
            
            ind2 = (JJ-1+1)*nn + KK;
            if ((ind2 > 0) && (ind2 <= length(vx)) && (mfi(ind2) > 0))
                A(mfi(ind1),mfi(ind2)) = 1;
            end
            ind2 = (JJ-1-1)*nn + KK;
            if ((ind2 > 0) && (ind2 <= length(vx)) && (mfi(ind2) > 0))
                A(mfi(ind1),mfi(ind2)) = 1;
            end
            ind2 = (JJ-1)*nn + KK+1;
            if ((ind2 > 0) && (ind2 <= length(vx)) && (mfi(ind2) > 0))
                A(mfi(ind1),mfi(ind2)) = 1;
            end
            ind2 = (JJ-1)*nn + KK-1;
            if ((ind2 > 0) && (ind2 <= length(vx)) && (mfi(ind2) > 0))
                A(mfi(ind1),mfi(ind2)) = 1;
            end
            
        end
        
    end
    
    f = sparse(A)\b;
    
    v(mif, n+1) = f;
    
    
    % normalize
    [~, vpos] = max(abs(real(v(:, 1)))); % max of zeroth mode
    %vpos = 1;
    v(:, n+1) = v(:, n+1) / v(vpos, n+1);
    
end
vx = vx(mfi>0, :);
v = v(mfi>0, :);

%scatter3(vx(:, 1), vx(:, 2), v(:, 15), 'r.')
if (0)
sgn = sign(vx(:, 1));
sgn(sgn==0) = 1;
cr = sgn.*(vx(:, 1).^2 + vx(:, 2).^2).^0.5;
cr = cr / 0.23;
plot(cr, real(v(:, 16)), '.')
hold on
plot(cr, imag(v(:, 16)), '.')
end

% complex amplitudes
b = load('/home/gholami/work/tests/fig18.csv'); % max velocity over time
U = compute_complex_amplitudes(b, 1000);

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