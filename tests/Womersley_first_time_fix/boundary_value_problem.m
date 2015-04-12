clear
clc

msh2d = 0;

if (msh2d)

r = 0.23;
c = 0;
dx = 0.04*.23;
na = floor(2*pi*r/dx) + 1;
da = 2*pi / na;

for k = 1:na
    c = c + 1;
    
    a = da * (k-1);
    x(c,1) = r * cos(a);
    y(c,1) = r * sin(a);
    
end

[vx,fvlist] = mesh2d([x y]);

else


r = 0.23;
c = 1;
x(1,1) = 0;
y(1,1) = 0;
dx = 0.09*.23;

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

end


%%% full connvectivity lists
% compute max degree
maxdeg = 0;
for k = 1:length(vx)
    maxdeg = max(maxdeg, numel(fvlist) - length(fvlist(abs(fvlist-k)>0)));
end

vflist = zeros(length(vx), maxdeg);
for k = 1:length(fvlist)
    for j = 1:3
        vflist(fvlist(k,j), :) =[k , vflist(fvlist(k,j),1:end-1)];
    end
end

vvlist = zeros(length(vx), maxdeg);
for k = 1:length(vx)
    temp = zeros(1,maxdeg);
    c = 1;
    while ((c <= maxdeg) && (vflist(k,c) > 0))
        for j = 1:3
            vi = fvlist(vflist(k,c), j);
            if ((vi ~= k) && (min(abs(temp - vi)) > 0)) % if not already counted
                temp = [vi, temp(1:end-1)];
            end
        end
        c = c + 1;
    end
    vvlist(k, :) = temp;
end

if (msh2d)
%%% boundary vertices
cc = 0;
for k = 1:length(vx)
    
    c = 1;
    while((c <= maxdeg) && (vvlist(k, c) > 0))
        vi = vvlist(k, c);
        
        common_f = vflist(k, ismember(vflist(k, :), vflist(vi, :)));
        common_f = common_f(common_f>0);
        
        if (length(common_f) == 1)
            cc = cc + 2;
            fe(cc-1:cc) = [k, vi];
        end
        
        c = c + 1;
    end
    
end
fe = unique(fe);
plot(vx(fe, 1), vx(fe, 2), 'ro')

end

%%% mapping indices between vx(full vertices) and vint(only interior)
c = 0;
mfi = -1 * ones(length(vx), 1);
for k = 1:length(vx)
    
    if (min(abs(fe - k)) > 0) % if NOT boundary
        
        c = c + 1;
        mfi(k) = c; % map full to interior
        mif(c) = k; % map interior to full
    end
    
end


%%% compute vertex mass and weights for Laplacian

% weight
ag = zeros(size(fvlist));
for k = 1:length(fvlist)
    for j = 1:3
        v1 = vx(fvlist(k, mod(j,3)+1), :) - vx(fvlist(k, j), :);
        v2 = vx(fvlist(k, mod(j+1,3)+1), :) - vx(fvlist(k, j), :);
        
        ag(k,j) = acos(dot(v1, v2) / norm(v1) / norm(v2));
        
    end
end

% mass
mass = zeros(length(vx), 1);
for k = 1:length(vx)
    
    if (min(abs(fe - k)) > 0) % i.e. if not a boundary vertex
        
        temp = zeros(maxdeg, 2);
        count = 0;
        c = 1;
        
        while((c <= maxdeg) && (vflist(k, c) > 0))
            count = count + 1;
            temp(count, 1) = mean(vx(fvlist(vflist(k, c), 1:3), 1));
            temp(count, 2) = mean(vx(fvlist(vflist(k, c), 1:3), 2));
            
            c = c + 1;
        end
        
        mass(k) = polyarea(temp(1:count, 1), temp(1:count, 2));
        
    end
    
end

%%% solve the boundary value problem
mu = 0.04;
T  = 0.92;
N = 24;
v = zeros(length(vx), N+1);
v(fe, :) = 0;
tic
for n = 0:N
    
    wn = 2*pi*n/T;
    kn2 = -1i * wn / mu;
    %kn2 = kn2 * 0.4150^2 *1.5;
    kn2 = kn2 * 4.5;
    %fold = zeros(length(vx)-length(fe), 1);
    %res = 1;
    A = zeros(length(vx)-length(fe), length(vx)-length(fe));
    b = zeros(length(vx)-length(fe), 1);
    
    %while (res > 0.1)
    
    
    for ind1 = 1:length(vx)
        
        if (min(abs(fe - ind1)) > 0) % i.e. if not a boundary vertex
            
            
            c = 1;
            while((c <= maxdeg) && (vvlist(ind1, c) > 0))
                ind2 = vvlist(ind1, c);
                
                
                common_f = vflist(ind1, ismember(vflist(ind1, :), vflist(ind2, :)));
                common_f = common_f(common_f>0);
              
                
                alpha1 = 0;
                for jj = 1:3
                    if ((fvlist(common_f(1), jj) ~= ind1) && (fvlist(common_f(1), jj) ~= ind2))
                        alpha1 = ag(common_f(1), jj);
                    end
                end
                alpha2 = 0;
                for jj = 1:3
                    if ((fvlist(common_f(2), jj) ~= ind1) && (fvlist(common_f(2), jj) ~= ind2))
                        alpha2 = ag(common_f(2), jj);
                    end
                end
                
                w = -(cot(alpha1) + cot(alpha2)) / 2;
                
                if (mfi(ind2) ~= -1)
                    A(mfi(ind1),mfi(ind2)) = -w;
                end % boundary is zero, so no else is required
                
                A(mfi(ind1),mfi(ind1)) = A(mfi(ind1),mfi(ind1)) + w;
                
                c = c + 1;
                
            end
            
            
            A(mfi(ind1),mfi(ind1)) = A(mfi(ind1),mfi(ind1)) + mass(ind1) * (kn2);
            b(mfi(ind1)) = -mass(ind1);
            
        end
        
    end
    
    f = sparse(A)\b;
    
    v(mif, n+1) = f;
    
    
    % normalize
    [~, vpos] = max(abs(real(v(:, 1)))); % max of zeroth mode
    %vpos = 1;
    v(:, n+1) = v(:, n+1) / v(vpos, n+1);
    
    %rv = real(v(:, n+1));
    %iv = imag(v(:, n+1));
    %v(:, n+1) = rv .* exp(-1i*iv);
    
    
end
toc
if (0)
cr = sign(vx(:, 1)).*(vx(:, 1).^2 + vx(:, 2).^2).^0.5;
cr = cr / 0.23;
ppp = polyfit(cr,v(:, 2),9);
cr = sort(cr);
pvv = polyval(ppp, cr);
plot(cr, real(pvv),'linewidth', 2);
hold on 
plot(cr, imag(pvv),'linewidth', 2);
end

someu = ('/home/gholami/work/tests/sazonov_tableCI');
U = someu(:, 3) .* exp(1i*someu(:, 4));

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
    u = u + 1 * real(U(n+1) * v(:, n+1) * exp(1i * wn * t(it)));
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