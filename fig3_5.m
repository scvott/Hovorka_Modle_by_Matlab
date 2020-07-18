clc
clear
close all
a=0;
b=1500;
step = 0.1;
t=a:step:b;
%parameter for Meals
tau_d = 40;
Ag = 0.8;%CHO to glucose utilization
Dm1(1) = 0;
Dm2(1) = 0;
Mwg = 180;%Molecular weight of glouse

BW = 70;%weight of the patient
u(1) = 12.9127;% ubasal is related to BW
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
uset = 0;

Dm1_2(1) = 0;
Dm2_2(1) = 0;
Q1_2(1) = 0.8*BW;
Q2_2(1) = -0.2292*BW+4.5307*u(1);
S1_2(1) = tau_s*u(1);
S2_2(1) = tau_s*u(1);
I_2(1) 	= u(1)/(0.01656*BW);
x1_2(1) = x1(1);
x2_2(1) = x2(1);
x3_2(1) = x3(1);
u_2(1)  = u(1);
uset_2 = 600;

Dm1_3(1) = 0;
Dm2_3(1) = 0;
Q1_3(1) = 0.8*BW;
Q2_3(1) = -0.2292*BW+4.5307*u(1);
S1_3(1) = tau_s*u(1);
S2_3(1) = tau_s*u(1);
I_3(1) 	= u(1)/(0.01656*BW);
x1_3(1) = x1(1);
x2_3(1) = x2(1);
x3_3(1) = x3(1);
u_3(1)  = u(1);
uset_3 = 450;
for i = 1 : length(t)
    [d_cho(i), Ug(i), Dm1(i+1), Dm2(i+1),G(i), Q1(i+1), Q2(i+1), Ui(i), S1(i+1), S2(i+1), I(i+1), x1(i+1), x2(i+1), x3(i+1), u(i+1)] = Hovorka_uset(step, tau_d, Ag, Mwg, F01, EGP, Vg, k12, Vi, tau_s, ka1, ka2, ka3, ke, Sit, Sid, Sie, Dm1(i), Dm2(i), Q1(i), Q2(i), S1(i), S2(i), I(i), x1(i), x2(i), x3(i), u(i), uset, t(i));
	
	[d_cho_2(i), Ug_2(i), Dm1_2(i+1), Dm2_2(i+1),G_2(i), Q1_2(i+1), Q2_2(i+1), Ui_2(i), S1_2(i+1), S2_2(i+1), I_2(i+1), x1_2(i+1), x2_2(i+1), x3_2(i+1), u_2(i+1)] = Hovorka_uset(step, tau_d, Ag, Mwg, F01, EGP, Vg, k12, Vi, tau_s, ka1, ka2, ka3, ke, Sit, Sid, Sie, Dm1_2(i), Dm2_2(i), Q1_2(i), Q2_2(i), S1_2(i), S2_2(i), I_2(i), x1_2(i), x2_2(i), x3_2(i), u_2(i), uset_2, t(i));
	
	[d_cho_3(i), Ug_3(i), Dm1_3(i+1), Dm2_3(i+1),G_3(i), Q1_3(i+1), Q2_3(i+1), Ui_3(i), S1_3(i+1), S2_3(i+1), I_3(i+1), x1_3(i+1), x2_3(i+1), x3_3(i+1), u_3(i+1)] = Hovorka_uset(step, tau_d, Ag, Mwg, F01, EGP, Vg, k12, Vi, tau_s, ka1, ka2, ka3, ke, Sit, Sid, Sie, Dm1_3(i), Dm2_3(i), Q1_3(i), Q2_3(i), S1_3(i), S2_3(i), I_3(i), x1_3(i), x2_3(i), x3_3(i), u_3(i), uset_3, t(i));
	
	time(i) = t(i)/60; %Unit of time is hour
    gbasul(i) = 5;     %g_based
end

figure(1)
hold on;
xlabel('Time (h)'); 
ylabel('D (g CHO/min)');
plot(time,d_cho(1:length(t)));
plot(time,d_cho_2(1:length(t)),'r--');
hold off;
legend('Meal_1', 'Meal_2'); 


figure(2)
hold on;
xlabel('Time (h)'); 
ylabel('G (mmol/min)');
plot(time,G(1:length(t)));
plot(time,G_2(1:length(t)),'r--');
plot(time,G_3(1:length(t)),'g--');
plot(time,gbasul(1:length(t)),'m');
hold off;
legend('Basal insulin', 'Bolus = 600', 'Bolus = 450', 'Normal glucose'); 

figure(3)
hold on;
xlabel('Time (h)');
ylabel('Ug (mmol/min)');
plot(time,Ug(1:length(t)));
plot(time,Ug(1:length(t)),'r--');
hold off;
legend('Glucose absorption_1', 'Glucose absorption_2');