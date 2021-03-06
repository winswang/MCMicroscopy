clear all;close all;
%% 1
x = -3:0.001:3;
a = [0.01,0.1,1,10];
figure(1);
for i = 1:length(a)
    y = sinc(x/a(i))/a(i);
    plot(x,y);hold on
end
xlabel('x');ylabel('sinc(x/a)/a');
legend('a = 10','a = 1', 'a = 0.1', 'a = 0.01');

%% 3
omega_x = 1;
omega_y = 1;
t = 0:pi/100:5*pi;
% case 1: Ex0 = Ey0
Ex0 = 1;
Ey0 = Ex0;
% case 1.1 delta = 0
delta = 0;
Ex = Ex0*sin(omega_x*t);
Ey = Ey0*sin(omega_y*t + delta);
figure(2);
subplot(231);
plot3(Ex,t,Ey);grid on
title('\delta = 0, E_{x0} = E_{y0}, linear');
xlabel('E_x');ylabel('t');zlabel('E_y');
axis([-1 1 0 16 -1 1]);
% case 1.2 delta = pi/2
delta = pi/2;
Ex = Ex0*sin(omega_x*t);
Ey = Ey0*sin(omega_y*t + delta);
subplot(232);
plot3(Ex,t,Ey);grid on
title('\delta = \pi/2, E_{x0} = E_{y0}, circular');
xlabel('E_x');ylabel('t');zlabel('E_y');
axis([-1 1 0 16 -1 1]);
% case 1.3 delta = pi/3
delta = pi/3;
Ex = Ex0*sin(omega_x*t);
Ey = Ey0*sin(omega_y*t + delta);
subplot(233);
plot3(Ex,t,Ey);grid on
title('\delta = \pi/3, E_{x0} = E_{y0}, circular');
xlabel('E_x');ylabel('t');zlabel('E_y');
axis([-1 1 0 16 -1 1]);
% case 2: Ex0 not eq Ey0
Ey0 = 0.5*Ex0;
% case 2.1 delta = 0
delta = 0;
Ex = Ex0*sin(omega_x*t);
Ey = Ey0*sin(omega_y*t + delta);
subplot(234);
plot3(Ex,t,Ey);grid on
title('\delta = 0, E_{x0} \neq E_{y0}, linear');
xlabel('E_x');ylabel('t');zlabel('E_y');
axis([-1 1 0 16 -1 1]);
% case 2.2 delta = pi/2
delta = pi/2;
Ex = Ex0*sin(omega_x*t);
Ey = Ey0*sin(omega_y*t + delta);
subplot(235);
plot3(Ex,t,Ey);grid on
title('\delta = \pi/2, E_{x0} \neq E_{y0}, elliptical');
xlabel('E_x');ylabel('t');zlabel('E_y');
axis([-1 1 0 16 -1 1]);
% case 1.3 delta = pi/3
delta = pi/3;
Ex = Ex0*sin(omega_x*t);
Ey = Ey0*sin(omega_y*t + delta);
subplot(236);
plot3(Ex,t,Ey);grid on
title('\delta = \pi/3, E_{x0} \neq E_{y0}, elliptical');
xlabel('E_x');ylabel('t');zlabel('E_y');
axis([-1 1 0 16 -1 1]);

%% 4
figure(3);
theta = 0:0.01:pi;
g = [1,0.8,0.2,0.01];
for i = 1:length(g)
    gi = g(i);
    p = (1-gi^2)./(4*pi*(1+gi^2-2*gi*cos(theta)).^(3/2));
    subplot(2,2,i);
    plot(theta,p);xlabel('\theta');ylabel('p(\theta)');
    title_label = strcat('g=',num2str(gi));
    title(title_label);
    axis([0 pi -0.1 4]);
end

%% 5
figure(4);
theta = 0:0.01:2*pi;
p = 3/4*(1+cos(theta).^2);
polar(theta,p);