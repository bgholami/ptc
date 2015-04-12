function Y = compute_complex_amplitudes(b, N)

x = linspace(b(1,1), b(end,1), N+1)';
y = interp1(b(:, 1), b(:, 2), x);

Nk = N;
Y = zeros(Nk+1, 1);
for k = 0:Nk
   
    for n = 0:N-1
        
        Y(k+1, 1) = Y(k+1, 1) + y(n+1) * exp(-1i*2*pi*k*n/N) / N;
        
    end
    
end

end