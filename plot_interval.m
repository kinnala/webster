function h = plot_interval(x1,x2,y1,y2,Dy1,Dy2, delta)

a = -(-(Dy1*x1) - Dy2*x1 + Dy1*x2 + Dy2*x2 + 2*y1 - 2*y2)/(x1-x2).^3;
b = -((Dy1*x1.^2 + 2*Dy2*x1.^2 + Dy1*x1*x2 - Dy2*x1*x2 - 2*Dy1*x2.^2 - Dy2*x2.^2 - 3*x1*y1 - 3*x2*y1 + 3*x1*y2 + 3*x2*y2)/(x1 - x2).^3);
c = -(-(Dy2*x1) + Dy1*x2)/(x1 - x2) - (3*x1*x2*(-(Dy1*x1) - Dy2*x1 + Dy1*x2 + Dy2*x2 + 2*y1 - 2*y2))/(x1 - x2).^3;
d = -((Dy2*x1.^3*x2 + Dy1*x1.^2*x2.^2 - Dy2*x1.^2*x2.^2 - Dy1*x1*x2.^3 - 3*x1*x2.^2*y1 + x2.^3*y1 - x1.^3*y2 + 3*x1.^2*x2*y2)/(x1 - x2).^3);

xs = x1:delta:x2;
h = plot(xs, a*xs.^3 + b*xs.^2 + c*xs + d);

end
