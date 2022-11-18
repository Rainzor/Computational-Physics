h =@(x,y) -2*(x.^2+y.^2)+1/2*(x.^4+y.^4)+1/2*(x-y).^4;
F = @(x,y,b) exp(-b.*h(x,y));
G = @(x,y,u,s) exp(-((x-u).^2+(y-u).^2)/(2*s^2))/(2*pi*s^2);

beta=5;
sigma = 1/(sqrt(2*beta));
% 
% A=importdata("data_beta_5.csv");
% d = A.data;
% scatplot(d(:,1),d(:,2),'circles', sqrt((range(d(:, 1))/30)^2 + (range(d(:,2))/30)^2), 100, 5, 1, 8);
% title(['\beta = ' num2str(beta) '的抽样图像'])

figure
[x,y] = meshgrid(-2:0.05:2);
z = F(x,y,beta);
surf(x,y,z)
set(gca,'ygrid','on','gridlinestyle','--','Gridalpha',0.4)
view(0,90)
title(['\beta = ' num2str(beta) '的分布图像'])
xlabel("x")
ylabel("y")
zlabel("F(x,y)")

% figure
% [x,y] = meshgrid(-6:0.05:6);
% z = G(x,y,sqrt(2),sigma)/2+G(x,y,-sqrt(2),sigma)/2;
% mesh(x,y,z)
% view(45,30)
% title(['\sigma = ' num2str(sigma) '的建议分布图像'])
% xlabel("x")
% ylabel("y")
% zlabel("T(x,y)")
% 
colormap(parula)
colorbar