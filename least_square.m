clc;
clear all;
close all;

x = meshgrid(0:50);
z = meshgrid(0:50);

y=0.5*z+0.4*x+25;
% y=0.4*x+25;
% noise=10*max(y)*rand(1,length(x),1)/100;
% y=y+noise;
G = [ones(length(x),1) x' z'];
% G = [ones(length(x),1) x' (x.*x)'];

m = inv(G'*G)*G'*y';
d1 = G*m;
figure(1)
hold on
plot3 (x,z,y,':bo','MarkerFaceColor','g','LineWidth',2);
mesh (x,z,d1);
legend ('Data','Regresi Linear',2)

% outlier
y2 = y;
y2(1) = y(1)*5*abs(y(1)-max(y))/max(y);
m2 = inv(G'*G)*G'*y2';
d2 = G*m2;
figure(2)
hold on
plot3 (x,z,y2,':bo','MarkerFaceColor','g','LineWidth',2);
mesh (x,z,d2);
legend ('Data dengan pencilan','Regresi Linear',2)

% Weighted Least Square
W = eye (length(x));
W(1,1) = 0.1;
m3 = inv(G'*W*G)*G'*W*y2';
d3 = G*m3;
figure(3)
hold on
plot3 (x,z,y2,':bo','MarkerFaceColor','g','LineWidth',2);
mesh (x,z,d3);
legend ('Data dengan pencilan','Regresi Linear',2)

%comparison
% figure (4)
% hold on
% plot (x,z,y,':bo','MarkerFaceColor','g','LineWidth',2);
% plot (x,z,d1,'b','LineWidth',2);
% plot (x,z,d2,'k','LineWidth',2);
% plot (x,z,d3,'r','LineWidth',2);
% legend ('Data','Regresi','Regresi(Outlier)','Regresi berbobot',2)
% 
