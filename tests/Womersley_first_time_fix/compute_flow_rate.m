function flowrate = compute_flow_rate(x, y)


s = 0;
for i = 1:length(x)-1
    tdx = x(i+1) - x(i);
    ty = (y(i+1) + y(i) ) / 2;
    s = s + tdx * ty;
end

flowrate = s;