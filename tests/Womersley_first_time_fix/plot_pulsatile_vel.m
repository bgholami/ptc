clear
clc


fol = '/scratch/run/testtaper3dinout/work/';
fol = '/home/gholami/work/port/smBCtest/work/';
name = 'mcf_particles00';
extension = 'out.crop';
num_digits = 6;
start_step = 30000;
end_step   = 435000;
jump = 1000;

dt =  0.10240775E-03;

% prepare for reading files one by one
num_steps = (end_step - start_step) / jump + 1;

step = start_step;

for i=1:num_steps
    
    step_l = length(num2str(step));
    step_str = '';
    for j=1:num_digits-step_l
        step_str = [step_str '0'];
    end
    step_str = [step_str num2str(step)];
    
    filename = [fol name step_str '.' extension];
    
    file_exist = exist(filename, 'file');
    
    if (file_exist == 2)
        filename
        
        a = load(filename);
        r = (a(:,2).^2 + a(:,3).^2).^0.5;
        p2 = polyfit(r, a(:, 4), 2);
        p3 = polyfit(r, a(:, 4), 3);

        v(i) = (polyval(p2, 0) + polyval(p3, 0)) / 2;
        %v(i) = polyval(p2, 0);
        t(i) = step * dt;
        
        step = step + jump;
    else
        num_steps = i-1;
        break
    end
end

for i = 1:1
    v = smooth(v);
end
plot(t, v * (2*r) / (pi * r^2), '--r','LineWidth',2)
plot(t/9.1451, v* 65 / .5852, '--r','LineWidth',2) % Re
hold on 

if (1)
    b = load('Hinds_pulsatile.csv'); % flow rate (ml/sec)
    
    r = 0.635 / 2; % cm
    tt = b(:, 1); % sec
    ttd = mean(tt(2:end) - tt(1:end-1));
    ttl = tt(end) - tt(1);
    vmean = b(:, 2) / (pi * r^2);  % cm/sec (i.e. vmean = c / area)
    vmax = 3/2 * vmean; % poiseuille
    % or just
    %b = load('hinds_vmax.data');

    
%     [~, mina] = min(a);
%     [~, minb] = min(vmax);
%     tt = tt + t(mina) - tt(minb);

    [~, ti] = min(vmax);
    vmax = [vmax(ti:end) ; vmax(1:ti-1)];
    tt = [tt(ti:end) ; tt(1:ti-1)+ttl+ttd];
    tt = tt - tt(1) + 30000*dt;
    
    av = vmax;
    at = tt;
    for i = 1:5
    av = [av ; vmax];
    at = [at ; tt + i*ttl + ttd];
    end
    hold on
    plot(at/9.1451, av * 64 / 1.13, '-k','LineWidth',2) % Re
end