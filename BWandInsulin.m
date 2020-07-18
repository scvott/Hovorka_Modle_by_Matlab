clc
clear
close all
%u(1) = 12.9127;% ubasal is related to BW
syms k;
% "solve" can solve the equation in matlab
t = 40:1:120;
for i = 1:120
    BW(i) = i+40;
    u_s = solve(0.2242*k^(2)+0.0405*k*BW(i)-0.0151*BW(i)^2 == 0,k,'Real',true);
    u_list = double(u_s);
    u(i) = u_list(2);
end


figure(1)
axis([40, 120, 3, 32]);
hold on;
title("Basal insulin infusion rate as function of the body weight");
xlabel('Body Weight (kg)'); 
ylabel('Basal insulin (mU/msin)');
plot(u);
hold off;

