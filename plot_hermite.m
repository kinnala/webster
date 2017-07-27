function [ h ] = plot_hermite( u,p )
N=length(p);
for jtr=1:(N-1)
    hold on;
    h=plot_interval(p(jtr),p(jtr+1),u(jtr),u(jtr+1),u(jtr+N),u(jtr+1+N),(p(jtr+1)-p(jtr))/10);
end
end

