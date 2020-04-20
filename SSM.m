clear 
clc
close

load("Satellite_1.mat")
mu =  3.98604419e14;
%a_0 = mu /(sqrt(sum(x_0(1:3).^2))).^2;

dt = 1; % time setp
N = 100000;
x = zeros(6,N);
x(:,1) = x_0; %add initial conditions to state matrix
x_new = x_0';
% t = 0;
for i = 2 : N
%     x_current = x_new;
%     x_new = x_current + [x_current(4:6); -(x_current(1:3)) ./(sqrt(sum(x_current(1:3).^2))) .* (mu / (sqrt(sum(x_current(1:3).^2))).^2)]*dt;
      t = t + dt;
%     hold on
%     scatter3(x_new(1,:),x_new(2,:),x_new(3,:));
    x(:,i) = x(:,i-1) + [x(4:6,i-1); -(x(1:3,i-1)) ./(sqrt(sum(x(1:3,i-1).^2))) .* (mu / (sqrt(sum(x(1:3,i-1).^2))).^2)] *dt;
end
t_days = t/(3600*24);
scatter3(x(1,:),x(2,:),x(3,:))
hold on
