function [v, vx] = solve_boundary_value_problem(vx, a, N, T, mu, nn, dx)


%%% mapping indices between vx(full vertices) and vint(only interior)
c = 0;
mfi = -1 * ones(length(vx), 1);
for k = 1:length(vx)
     % this is for a circular problem, for an arbitrary boundary patch,
     % this mapping should be provided
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

end