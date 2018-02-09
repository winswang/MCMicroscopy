close all
for n = 1:3
    xStart = 0;
    yStart = 0;
    xPoint = zeros(1,20);
    yPoint = zeros(1,20);
    for m = 1:20
        x1 = rand(1);
        s = -100*log(x1);
        x2 = rand(1);
        alpha = 2*pi*x2;
        xPoint(m) = s*cos(alpha);
        yPoint(m) = s*sin(alpha);
    end
    x = [xStart,xPoint];
    y = [yStart,yPoint];
    figure;
    plot(x,y);
    xlabel('x');ylabel('y');axis square
end