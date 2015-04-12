clear
clc

a = load('/home/gholami/work/tests/sazonov_tableCI');
b = load('/home/gholami/Downloads/fig18.csv');

N = 24;

x = linspace(b(1,1), b(end,1), N)';
y = interp1(b(:, 1), b(:, 2), x);
figure
plot(x, y)
dx = x(2) - x(1);
T = x(end) - x(1);

u = a(:, 3) .* exp(1i*a(:, 4));

for k = 0:N-1
   
    U(k+1,1) = sum(u.*exp(1i*a(:,2)*2*pi*x(k+1)));
    
end

hold on 
plot(x, real(U), 'r')

% Y = zeros(24+1, 1);
% for k = 0:24
%    
%     %Y(n+1,1) = sum(y.*exp(-1i*2*pi/T*[0:N]'*x(k´n+1)));
%     wn = 2*pi*k/T;
%     for j = 1:length(y)
%         
%         Y(k+1, 1) = Y(k+1, 1) + y(j) * exp(-1i*wn*x(j)) * dx;
%         
%     end
%     
% end

Nk = 24; % ??
Y = zeros(Nk+1, 1);
for k = 0:Nk
   
    for n = 0:N-1
        
        Y(k+1, 1) = Y(k+1, 1) + y(n+1) * exp(-1i*2*pi*k*n/N) / N;
        
    end
    
end


figure
plot(abs(Y))
hold on 
plot(abs(u), 'r')

figure
plot(angle(Y))
hold on 
plot(angle(u), 'r')