K = @(a,b)(bsxfun(@min, a', b));
f = @(x)(exp(cosh((x + 2*x.^2+cos(x))./(3+sin(x.^3)))));

f = @(x)(exp(-sin(3*(6*x-3)).^2 - (6*x-3).^2));

x = linspace(0,1,500);
plot(x, f(x))

%%

figure(1); clf; hold on
title('Brownian motion Kernel', 'Fontsize', 25)

for i = 0:10

    X = ((0:2^i)/2^i)';
    Y = f(X);

    Q = K(X,X) + 10d-7 * eye(size(X,1));
    m = K(x,X') * (Q\Y); 
    sigma = 2*sqrt(diag(K(x,x)+ 10d-7 * eye(size(x,1)) - K(x,X') * (Q\K(X',x))));

    
    fill([x';flipud(x')],[m-sigma;flipud(m+sigma)],'b','linestyle','none');
    plot(x, m, 'k--'); 
    alpha(0.05)
    scatter(X,zeros(length(X),1), 'filled', 'red')
    axis([0 1 0 1])
    drawnow
    pause

end

hold off

%%

K = @(a,b) (0.5 * abs(bsxfun(@minus,a',b)) .* bsxfun(@min,a',b).^2 + (1/3) * bsxfun(@min,a',b).^3 );

figure(2); clf; hold on
title('Spline Kernel', 'Fontsize', 25)

for i = 0:10

    X = ((0:2^i)/2^i)';
    Y = f(X);

    Q = K(X,X) + 10d-7 * eye(size(X,1));
    m = K(x,X') * (Q\Y); 
    sigma = 2*sqrt(diag(K(x,x)+ 10d-7 * eye(size(x,1)) - K(x,X') * (Q\K(X',x))));

    
    fill([x';flipud(x')],[m-sigma;flipud(m+sigma)],'b','linestyle','none');
    plot(x, m, 'k--'); 
    alpha(0.25)
    scatter(X,zeros(length(X),1), 'filled', 'red')
    axis([0 1 0 1])
    drawnow
    pause

end