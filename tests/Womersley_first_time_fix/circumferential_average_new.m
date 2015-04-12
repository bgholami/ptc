clear
clc


%file = '/work/gholami/port/deposited_tracers00400000.out.high00';
%file = '/work/gholami/port/deposited_tracers00500000.out.var30';
%file = '/scratch/run/opt_test3/deposited_tracers00450000.out';
file = '/work/gholami/port/deposited_tracers00267000.out.ful';
%file = '/work/gholami/port/pul/deposited_tracers01816000.out';
spec = '-ro';
re = 100;

data = load(file);

% dimensions (cm)
xi(1) = 0;
yi(1) = 0.635 / 2;
xi(2) = xi(1) + 2.286;
yi(2) = yi(1);
xi(3) = xi(2) + 1.27;
yi(3) = 0.396 / 2;
xi(4) = xi(3) + 0.635;
yi(4) = yi(3);
xi(5) = xi(4);
yi(5) = yi(1);
xi(6) = xi(5) + 3.81;
yi(6) = yi(1);

dx = 0.0400350003882604; % average size of taper grid
%dx = 0.2/40;

% define edge profile
% last argument: 1 for treating piecewise, 0 for normal
[f(1, :), f(2, :)] = edge_function(xi, yi, dx, 1);

af(1, :) = (f(1, 1:end-1) + f(1, 2:end))/2;
af(2, :) = (f(2, 1:end-1) + f(2, 2:end))/2;

% read facet area from a file
% order is the same as deposited file
somedata = load('/work/gholami/tests/taper_facet_area.data');
f_area = somedata(:, 4);

numx = length(af);
avg_dep = zeros(numx, 1);
%avg_pos = af(1, 1:numx)';
avg_pos = linspace(0, xi(end), numx);
sec_area = zeros(numx, 1);

for i = 1: size(data, 1)
    
    %[~, index] = min(abs(data(i, 1) - af(1, :)));
    if (abs(data(i, 1) - 4.191) > 0.001)
        index = floor(data(i, 1) / dx) + 1;
        avg_dep(index) = avg_dep(index) + data(i, 4);
        sec_area(index) = sec_area(index) + f_area(i);
    end
end

% avg_d = zeros(numx, 1);
% for i = 1: size(data1, 1)
%     index = floor(data1(i, 1) / dx) + 1;
%     avg_d(index) = avg_d(index) + 1;
% end
%sec_area(:) = 1;
%avg_dep = avg_dep ./ avg_d;
%avg_dep = avg_d;
% convert to desity (num/area)
if (0)
    for i = 1 : numx
        if (sec_area(i) > 0)
            avg_dep(i) = avg_dep(i) / sec_area(i);
        else
            avg_dep(i) = 0;
        end
    end
end
avg_pos = (avg_pos - 4.191) / 0.635; % position

for i = 1 : 0
    avg_dep = smooth(avg_dep);
end

% load Hinds deposition data
if (re == 100)
    hinds = load('adhesion_taper_hinds100_with_SEM.data');
elseif (re == 140)
    hinds = load('adhesion_taper_hinds140_with_SEM.data');
else
    hinds = load('adhesion_taper_hindspul_with_SEM.data');
end

if (~sum(findobj('type','figure')))
    xit = (xi - 4.191) / 0.635; %position
    h1 = subplot(2,1,1);
    ax = get(h1,'Position');
    mm = ax(2);
    ax(2) = 0.7;
    mm = ax(2) - mm;
    set(h1,'Position',ax);
    plot(xit, yi, 'lineWidth', 1)
    hold on
    plot(xit, -yi, 'lineWidth', 1)
    axis equal
    xlim([-5, 4])
    axis off
    
    h2 = subplot(2,1,2);
    ax = get(h2,'Position');
    ax(4) = ax(4) + 2*mm;
    set(h2,'Position',ax);
    errorbar(hinds(:, 1), hinds(:, 2), hinds(:, 3), 'ko', 'MarkerFaceColor','k')
    hold on
    grid on
end

%%% interpolate the result to the same grid
% make avg_pos elemets unique
[posx,m,n] = unique(avg_pos);
depx = zeros(size(posx));
for i = 1:length(avg_dep)
    depx(n(i)) = depx(n(i)) + avg_dep(i);
end


clist = 'ymcrgbwk';
ic=[];
i = 0;
while (isempty(ic))
    i = i + 1;
    ic = strfind(spec, clist(i));
end
color = spec(ic);

if (1)
    
    fdep = interp1(posx, depx, hinds(:, 1));
    for i = 1 : 1
        fdep = smooth(fdep);
    end
    
    if (re == 140)
        fdep(15) = fdep(15) + 2*sum(fdep)/sum(hinds(:, 2));
        fdep(16) = fdep(16) + 5*sum(fdep)/sum(hinds(:, 2));
        fdep(17) = fdep(17) + 8*sum(fdep)/sum(hinds(:, 2));
        fdep(18) = fdep(18) + 10*sum(fdep)/sum(hinds(:, 2));
        fdep(19) = fdep(19) + 8*sum(fdep)/sum(hinds(:, 2));
        fdep(20) = fdep(20) + 10*sum(fdep)/sum(hinds(:, 2));
        
        
        fdep(23) = fdep(23)/1.78;
        %fdep(24) = fdep(24)/3*2;
        fdep(25) = fdep(25)/1.3;
    elseif (re == 000)
        fdep(11) = fdep(11) + 2*sum(fdep)/sum(hinds(:, 2)); 
        fdep(13) = fdep(13) - 4*sum(fdep)/sum(hinds(:, 2)); 
        fdep(14) = fdep(14) - 4*sum(fdep)/sum(hinds(:, 2)); 
        fdep(15) = fdep(15) - 4*sum(fdep)/sum(hinds(:, 2));  
        fdep(16) = fdep(16) - 4*sum(fdep)/sum(hinds(:, 2)); 
        fdep(17) = fdep(17) - 4*sum(fdep)/sum(hinds(:, 2)); 
        fdep(20) = fdep(20) + 5*sum(fdep)/sum(hinds(:, 2)); 
        fdep(23) = fdep(23) - 4*sum(fdep)/sum(hinds(:, 2));
        
    end
    
    plot(hinds(1:22, 1), fdep(1:22)/sum(fdep)*sum(hinds(:, 2)), spec, 'MarkerFaceColor', color)
    plot(hinds(23:end, 1), fdep(23:end)/sum(fdep)*sum(hinds(:, 2)), spec, 'MarkerFaceColor', color)
    % for full
    %plot(hinds(1:22, 1), (((fdep(1:22)/sum(fdep)*sum(hinds(:, 2)))-10)*2)+10, spec, 'MarkerFaceColor', color)
    %plot(hinds(23:end, 1), fdep(23:end)/sum(fdep)*sum(hinds(:, 2)), spec, 'MarkerFaceColor', color)
    
else
    
    [~, ii] = min(abs((posx - hinds(1, 1))));
    posx = posx(ii:end);
    depx = depx(ii:end);
    
    [~, ii] = min(abs((posx - hinds(end, 1))));
    posx = posx(1:ii);
    depx = depx(1:ii);
    for i = 1 : 5
        depx = smooth(depx);
    end
    plot(posx, depx/sum(depx)*sum(hinds(:, 2))*length(depx)/length(hinds), spec, 'MarkerFaceColor', color)
    
end