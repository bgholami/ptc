clear
clc

xx = -1:.05:1; yy = abs(sqrt(xx));
[x,y] = pol2cart(xx,yy);
k = convhull(x,y);
hold on
plot(x(k),y(k),'r-',x,y,'b+')
x_test = rand(10,1);
y_test = rand(10,1);
x = [x' ; x(1)]; y = [y' ; y(1)];
IN = inpolygon(x_test,y_test,x, y);%x(k),y(k))
plot(x_test(IN),y_test(IN),'g*',x_test(~IN),y_test(~IN),'y*')
hold off