function flowrate = compute_flow_rate3D(x, y)
% assumes rotation around 0, i.e. r = abs(x)

s = 0;
for i = 1:length(x)-1
    dA = pi/2 * abs(x(i+1)^2 - x(i)^2);
    ty = (y(i+1) + y(i) ) / 2;
    s = s + dA * ty;
end

flowrate = s;