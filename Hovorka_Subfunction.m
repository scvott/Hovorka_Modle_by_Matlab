clc
close all
clear all
a=0;
b=1500;
step = 0.001;
t=a:step:b;
%parameter for Meals
tau_d = 40;
Ag = 0.8;%CHO to glucose utilization
Dm1(1) = 0;
Dm2(1) = 0;
Mwg = 180;%Molecular weight of glouse

BW = 70;%weight of the patient
%u(1) = 12.9127;% ubasal is related to BW
syms k;
% "solve" can solve the equation in matlab
u_s = solve(0.2242*k^(2)+0.0405*k*BW-0.0151*BW^(2) == 0,k,'Real',true);
u_list = double(u_s);
u(1) = u_list(2);
%parameter for Glucose
F01 = 0.00097*BW;
EGP = 0.0161*BW;
Vg = 0.16*BW;
k12 = 0.066;%transfer rate

%parameter for Insulin
Vi = 0.12*BW;
tau_s = 55;
x1(1) = 0.30898*u(1)/BW;
x2(1) = 0.04951*u(1)/BW;
x3(1) = 3.2206*u(1)/BW;
Q1(1) = 0.8*BW;
Q2(1) = -0.2292*BW+4.5307*u(1);
I(1) = u(1)/(0.01656*BW);

%parameter for Insulin Action
ka1 = 0.006;
ka2 = 0.06;
ka3 = 0.03;
ke = 0.138;
Sit = 51.2*10^(-4);
Sid = 8.2*10^(-4);
Sie = 520*10^(-4);	
kb1 = Sit*ka1;
kb2 = Sid*ka2;
kb3 = Sie*ka3;
S1(1) = tau_s*u(1);
S2(1) = tau_s*u(1);

eat = 10; %eatting 
for i = 1 : length(t)
    [d_cho(i), Ug(i), Dm1(i+1), Dm2(i+1),G(i), Q1(i+1), Q2(i+1), Ui(i), S1(i+1), S2(i+1), I(i+1), x1(i+1), x2(i+1), x3(i+1), u(i+1)] = Hovorka_sub(step, tau_d, Ag, Mwg, F01, EGP, Vg, k12, Vi, tau_s, ka1, ka2, ka3, kb1, kb2, kb3, ke, Sit, Sid, Sie, Dm1(i), Dm2(i), Q1(i), Q2(i), S1(i), S2(i), I(i), x1(i), x2(i), x3(i), u(i), eat, t(i));
	
	time(i) = t(i)/60; %Unit of time is hour
    gbasul(i) = 5;     %g_based
    g_high(i) = 8;     
    g_low(i) = 4;     
    
end

figure(1)
hold on;
xlabel('Time (h)'); 
ylabel('D (g CHO/min)');
plot(time,d_cho(1:length(t)));
hold off;
legend('Meal'); 

figure(2)
hold on;
xlabel('Time (h)'); 
ylabel('Insulin infusion rate (mU/min)');
plot(time,u(1:length(t)));
hold off;
legend('Insulin infusion rate');


figure(3)
hold on;
xlabel('Time (h)'); 
ylabel('G (mmol/min)');
plot(time,G(1:length(t)));
plot(time,gbasul(1:length(t)),'r--');
plot(time,g_high(1:length(t)),'m--');
plot(time,g_low(1:length(t)),'g--');
hold off;
legend('Glucose concentration', 'Setpoint', 'High blood sugar', 'Low blood sugar'); 

figure(4)
hold on;
xlabel('Time (h)'); 
ylabel('Ug (mmol/min)');
plot(time,Ug(1:length(t)));
hold off;
legend('Glucose absorption');