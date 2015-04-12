clear
clc

a = load('/home/gholami/Downloads/fig18.csv');

N = 24;

x = linspace(a(1,1), a(end,1), N+1);
y = interp1(a(:, 1), a(:, 2), x);
dx = x(2) - x(1);
T = x(end) - x(1);

n = 15;

u = 0;
for k = 0:N-1
    
    u = u + y(k+1) * exp(-i * 2 * pi * n * k / N);
    
    %t = (k - 1/2) * dx;
    %wn = 2 * pi * n / T;
    %u = u + (y(k) + y(k+1))/2 * exp(-i * wn * t);

end
u = u * dx;

abs(u)