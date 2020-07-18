clc
close all
clear all
a=0;
b=180;
step = 0.01;
t=a:step:b;
%parameter for Meals
tau_d = 40;
Ag = 0.8;%CHO to glucose utilization
Dm1(1) = 0;
Dm2(1) = 0;
Mwg = 180;%Molecular weight of glouse
eat_time = 5;

%parameter for Glucose
BW = 70;%weight of the patient
F01 = 0.00097*BW; 
Fc_01 = 0;
EGP = 0.0161*BW;
Vg = 0.16*BW;
k12 = 0.066;%transfer rate

%parameter for Insulin
Vi = 0.12*BW;
tau_s = 55;
x1(1) = 0.0295;
x2(1) = 4.72*10^(-3);
x3(1) = 0.3072;
Q1(1) = 412.55;
Q2(1) = 151.09;
u(1) = 6.678;
I(1) = 5.761;

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
S1(1) = 0;
S2(1) = 0;

for j = 1 : length(t)
    %Meals RungeKutta
    if t(j)<=eat_time
        d_cho(j)=10;
    else
		d_cho(j)=0;
	end
    D_meal(j) = 1000*d_cho(j)/Mwg;
    
	rkm11 = Ag*D_meal(j)-Dm1(j)/tau_d;
    rkm21 = Dm1(j)/tau_d-Dm2(j)/tau_d;
    
    rkm12 = Ag*D_meal(j)-(Dm1(j)+step*rkm11/2)/tau_d;
    rkm22 = (Dm1(j)+step*rkm11/2)/tau_d-(Dm2(j)+step*rkm21/2)/tau_d;
      
    rkm13= Ag*D_meal(j)-(Dm1(j)+step*rkm12/2)/tau_d;
    rkm23= (Dm1(j)+step*rkm12/2)/tau_d-(Dm2(j)+step*rkm22/2)/tau_d;
      
    rkm14= Ag*D_meal(j)-(Dm1(j)+step*rkm13)/tau_d;
    rkm24= (Dm1(j)+step*rkm13)/tau_d-(Dm2(j)+step*rkm23)/tau_d;
    
	Dm1(j+1) = Dm1(j)+step*(rkm11+2*rkm12+2*rkm13+rkm14)/6;
	Dm2(j+1) = Dm2(j)+step*(rkm21+2*rkm22+2*rkm23+rkm24)/6;
    
	Ug(j) = Dm2(j)/tau_d;
    
    %Glucose RungeKutta
	G(j) = Q1(j)/Vg;
	if G(j)>=4.5
		Fc_01 = F01;
	else
		Fc_01 = F01*G(j)/4.5;
	end
	if G(j)>=9
		Fr(j) = 0.003*(G(j)-9)*Vg;
	else
		Fr(j) = 0;
	end
	
	rkg11 = Ug(j)-Fc_01-Fr(j)-x1(j)*Q1(j)+k12*Q2(j)+EGP*(1-x3(j));
	rkg21 = x1(j)*Q1(j)-(k12+x2(j))*Q2(j);

	rkg12 = Ug(j)-Fc_01-Fr(j)-x1(j)*(Q1(j)+step/2*rkg11)+k12*(Q2(j)+step/2*rkg21)+EGP*(1-x3(j));
	rkg22 = x1(j)*(Q1(j)+step/2*rkg11)-(k12+x2(j))*(Q2(j)+step/2*rkg21);
	
	rkg13 = Ug(j)-Fc_01-Fr(j)-x1(j)*(Q1(j)+step/2*rkg12)+k12*(Q2(j)+step/2*rkg22)+EGP*(1-x3(j));
	rkg23 = x1(j)*(Q1(j)+step/2*rkg12)-(k12+x2(j))*(Q2(j)+step/2*rkg22);
	
	rkg14 = Ug(j)-Fc_01-Fr(j)-x1(j)*(Q1(j)+step*rkg13)+k12*(Q2(j)+step*rkg23)+EGP*(1-x3(j));
	rkg24 = x1(j)*(Q1(j)+step*rkg13)-(k12+x2(j))*(Q2(j)+step*rkg23);
	
	Q1(j+1) = Q1(j)+step*(rkg11+2*rkg12+2*rkg13+rkg14)/6;
	Q2(j+1) = Q2(j)+step*(rkg21+2*rkg22+2*rkg23+rkg24)/6;
	
	%Insulin
	rki11 = u(j)-S1(j)/tau_s;
	rki21 = (S1(j)-S2(j))/tau_s;
	
	rki12 = u(j)-(S1(j)+step/2*rki11)/tau_s;
	rki22 = ((S1(j)+step/2*rki11)-(S2(j)+step/2*rki21))/tau_s;
	
	rki13 = u(j)-(S1(j)+step/2*rki12)/tau_s;
	rki23 = ((S1(j)+step/2*rki12)-(S2(j)+step/2*rki22))/tau_s;
	
	rki14 = u(j)-(S1(j)+step*rki13)/tau_s;
	rki24 = ((S1(j)+step*rki13)-(S2(j)+step*rki23))/tau_s;
	
	S1(j+1) = S1(j)+step*(rki11+2*rki12+2*rki13+rki14)/6;
	S2(j+1) = S2(j)+step*(rki21+2*rki22+2*rki23+rki24)/6;
	
	Ui(j) = S2(j)/tau_s;
	
	%Insulin Action
	rkia1 = Ui(j)/Vi-ke*I(j);
	rkia2 = Ui(j)/Vi-ke*(I(j)+step/2*rkia1);
	rkia3 = Ui(j)/Vi-ke*(I(j)+step/2*rkia2);
	rkia4 = Ui(j)/Vi-ke*(I(j)+step*rkia3);
	I(j+1) = I(j)+step*(rkia1+2*rkia2+2*rkia3+rkia4)/6;
	
	rkx11 = -ka1*x1(j)+kb1*I(j);
	rkx21 = -ka2*x2(j)+kb2*I(j);
	rkx31 = -ka3*x3(j)+kb3*I(j);
	
	rkx12 = -ka1*(x1(j)+step/2*rkx11)+kb1*I(j);
	rkx22 = -ka2*(x2(j)+step/2*rkx21)+kb2*I(j);
	rkx32 = -ka3*(x3(j)+step/2*rkx31)+kb3*I(j);
	
	rkx13 = -ka1*(x1(j)+step/2*rkx12)+kb1*I(j);
	rkx23 = -ka2*(x2(j)+step/2*rkx22)+kb2*I(j);
	rkx33 = -ka3*(x3(j)+step/2*rkx32)+kb3*I(j);
	
	rkx14 = -ka1*(x1(j)+step*rkx13)+kb1*I(j);
	rkx24 = -ka2*(x2(j)+step*rkx23)+kb2*I(j);
	rkx34 = -ka3*(x3(j)+step*rkx33)+kb3*I(j);
	x1(j+1) = x1(j)+step*(rkx11+2*rkx12+2*rkx13+rkx14)/6;
	x2(j+1) = x2(j)+step*(rkx21+2*rkx22+2*rkx23+rkx24)/6;
	x3(j+1) = x3(j)+step*(rkx31+2*rkx32+2*rkx33+rkx34)/6;
	
	u(j+1) = u(j);
end

figure(1)
hold on;
xlabel('Time (min)'); 
ylabel('mmol/min');
plot(t,Ug(1:length(t)));
plot(t,G(1:length(t)));
hold off;
legend('Glucose asorption rate','Glucose'); 

figure(2)
hold on;
xlabel('Time (min)'); 
ylabel('mmol/min');
plot(t,Ui(1:length(t)));
plot(t,I(1:length(t)));
hold off;
legend('Insulin asorption rate','Insulin');

figure(3)
hold on;
xlabel('Time (min)'); 
ylabel('mmol/min');
plot(t,Q1(1:length(t)));
plot(t,Q2(1:length(t)));
hold off;
legend('Q1','Q2'); 


figure(4)
hold on;
xlabel('Time (min)'); 
ylabel('mmol/min');
plot(t,S1(1:length(t)));
plot(t,S2(1:length(t)));
hold off;
legend('S1','S2');


figure(5)
hold on;
xlabel('Time (min)'); 
ylabel('mmol/min');
plot(t,x1(1:length(t)));
plot(t,x2(1:length(t)));
plot(t,x3(1:length(t)));
hold off;
legend('x1','x2', 'x3');

figure(6)
hold on;
xlabel('Time (min)'); 
ylabel('mmol/min');
plot(t,Dm1(1:length(t)));
plot(t,Dm2(1:length(t)));
hold off;
legend('D1','D2');

