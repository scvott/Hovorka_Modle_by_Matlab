clc
clear
close all
a=0;
b=4500;
step = 0.001;
t=a:step:b;
day2 = 25*60;
day3 = 25*60*2;

%parameter for Meals
tau_d = 40;
Ag = 0.8;%CHO to glucose utilization
Dm1(1) = 0;
Dm2(1) = 0;
Mwg = 180;%Molecular weight of glouse

%0.2242*k^(2)+0.0405*k*BW-0.0151*Bw^(2)
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

for i = 1 : length(t)
    %Meals RungeKutta
    if t(i)>=210 && t(i)<=220
        d_cho(i)=12.42;
    elseif t(i)>=390 && t(i)<=410
        d_cho(i)=1.83;
    elseif t(i)>=600 && t(i)<=610
        d_cho(i)=3.45;
    elseif t(i)>=765 && t(i)<=785
        d_cho(i)=3.85;
    elseif t(i)>=900 && t(i)<=910
        d_cho(i)=3.45;
    elseif t(i)>=(210+day2) && t(i)<=(220+day2)
        d_cho(i)=12.42;
    elseif t(i)>=(390+day2) && t(i)<=(410+day2)
        d_cho(i)=1.83;
    elseif t(i)>=(600+day2) && t(i)<=(610+day2)
        d_cho(i)=3.45;
    elseif t(i)>=(765+day2) && t(i)<=(785+day2)
        d_cho(i)=3.85;
    elseif t(i)>=(900+day2) && t(i)<=(910+day2)
        d_cho(i)=3.45;
    elseif t(i)>=(210+day3) && t(i)<=(220+day3)
        d_cho(i)=12.42;
    elseif t(i)>=(390+day3) && t(i)<=(410+day3)
        d_cho(i)=1.83;
    elseif t(i)>=(600+day3) && t(i)<=(610+day3)
        d_cho(i)=3.45;
    elseif t(i)>=(765+day3) && t(i)<=(785+day3)
        d_cho(i)=3.85;
    elseif t(i)>=(900+day3) && t(i)<=(910+day3)
        d_cho(i)=3.45;
    else
		d_cho(i)=0;
    end
    % insullin inject
    if t(i)>=210 && t(i)<=215
         u(i) = 700;
    elseif t(i)>=390 && t(i)<=395
        u(i) = 250;
    elseif t(i)>=600 && t(i)<=605
        u(i) = 200;
    elseif t(i)>=765 && t(i)<=770
        u(i) = 490;
    elseif t(i)>=900 && t(i)<=905
        u(i) = 210;
    elseif t(i)>=(210+day2) && t(i)<=(215+day2)
        u(i) = 700;
    elseif t(i)>=(390+day2) && t(i)<=(395+day2)
        u(i) = 250;
    elseif t(i)>=(600+day2) && t(i)<=(605+day2)
        u(i) = 200;
    elseif t(i)>=(765+day2) && t(i)<=(770+day2)
        u(i) = 490;
    elseif t(i)>=(900+day2) && t(i)<=(905+day2)
        u(i) = 210;
    elseif t(i)>=(210+day3) && t(i)<=(215+day3)
        u(i) = 700;
    elseif t(i)>=(390+day3) && t(i)<=(395+day3)
        u(i) = 250;
    elseif t(i)>=(600+day3) && t(i)<=(605+day3)
        u(i) = 200;
    elseif t(i)>=(765+day3) && t(i)<=(770+day3)
        u(i) = 490;
    elseif t(i)>=(900+day3) && t(i)<=(905+day3)
        u(i) = 210;
    else
		u(i) = u(1);
    end
    
    D_meal(i) = 1000*d_cho(i)/Mwg;
    
	rkm11 = step*(Ag*D_meal(i)-Dm1(i)/tau_d);
    rkm21 = step*(Dm1(i)/tau_d-Dm2(i)/tau_d);
    
    rkm12 = step*(Ag*D_meal(i)-(Dm1(i)+rkm11/2)/tau_d);
    rkm22 = step*((Dm1(i)+rkm11/2)/tau_d-(Dm2(i)+rkm21/2)/tau_d);
      
    rkm13= step*(Ag*D_meal(i)-(Dm1(i)+rkm12/2)/tau_d);
    rkm23= step*((Dm1(i)+rkm12/2)/tau_d-(Dm2(i)+rkm22/2)/tau_d);
      
    rkm14= step*(Ag*D_meal(i)-(Dm1(i)+rkm13)/tau_d);
    rkm24= step*((Dm1(i)+rkm13)/tau_d-(Dm2(i)+rkm23)/tau_d);
    
	Dm1(i+1) = Dm1(i)+(rkm11+2*rkm12+2*rkm13+rkm14)/6;
	Dm2(i+1) = Dm2(i)+(rkm21+2*rkm22+2*rkm23+rkm24)/6;
    
	Ug(i) = Dm2(i)/tau_d;
    
    %Glucose RungeKutta
	G(i) = Q1(i)/Vg;
	if G(i)>=4.5
		Fc_01 = F01;
	else
		Fc_01 = F01*G(i)/4.5;
	end
	if G(i)>=9
		Fr(i) = 0.003*(G(i)-9)*Vg;
	else
		Fr(i) = 0;
	end
	
	rkg11 = step*(Ug(i)-Fc_01-Fr(i)-x1(i)*Q1(i)+k12*Q2(i)+EGP*(1-x3(i)));
	rkg21 = step*(x1(i)*Q1(i)-(k12+x2(i))*Q2(i));

	rkg12 = step*(Ug(i)-Fc_01-Fr(i)-x1(i)*(Q1(i)+rkg11/2)+k12*(Q2(i)+rkg21/2)+EGP*(1-x3(i)));
	rkg22 = step*(x1(i)*(Q1(i)+rkg11/2)-(k12+x2(i))*(Q2(i)+rkg21/2));
	
	rkg13 = step*(Ug(i)-Fc_01-Fr(i)-x1(i)*(Q1(i)+rkg12/2)+k12*(Q2(i)+rkg22/2)+EGP*(1-x3(i)));
	rkg23 = step*(x1(i)*(Q1(i)+rkg12/2)-(k12+x2(i))*(Q2(i)+rkg22/2));
	
	rkg14 = step*(Ug(i)-Fc_01-Fr(i)-x1(i)*(Q1(i)+rkg13)+k12*(Q2(i)+rkg23)+EGP*(1-x3(i)));
	rkg24 = step*(x1(i)*(Q1(i)+rkg13)-(k12+x2(i))*(Q2(i)+rkg23));
	
	Q1(i+1) = Q1(i)+(rkg11+2*rkg12+2*rkg13+rkg14)/6;
	Q2(i+1) = Q2(i)+(rkg21+2*rkg22+2*rkg23+rkg24)/6;
	
	%Insulin
	rki11 = step*(u(i)-S1(i)/tau_s);
	rki21 = step*((S1(i)-S2(i))/tau_s);
	
	rki12 = step*(u(i)-(S1(i)+rki11/2)/tau_s);
	rki22 = step*(((S1(i)+rki11/2)-(S2(i)+rki21/2))/tau_s);
	
	rki13 = step*(u(i)-(S1(i)+rki12/2)/tau_s);
	rki23 = step*(((S1(i)+rki12/2)-(S2(i)+rki22/2))/tau_s);
	
	rki14 = step*(u(i)-(S1(i)+rki13)/tau_s);
	rki24 = step*(((S1(i)+rki13)-(S2(i)+rki23))/tau_s);
	
	S1(i+1) = S1(i)+(rki11+2*rki12+2*rki13+rki14)/6;
	S2(i+1) = S2(i)+(rki21+2*rki22+2*rki23+rki24)/6;
	
	Ui(i) = S2(i)/tau_s;
	
	%Insulin Action
	rkia1 = step*(Ui(i)/Vi-ke*I(i));
	rkia2 = step*(Ui(i)/Vi-ke*(I(i)+rkia1/2));
	rkia3 = step*(Ui(i)/Vi-ke*(I(i)+rkia2/2));
	rkia4 = step*(Ui(i)/Vi-ke*(I(i)+rkia3));
	I(i+1) = I(i)+(rkia1+2*rkia2+2*rkia3+rkia4)/6;
	
	rkx11 = step*(-ka1*x1(i)+kb1*I(i));
	rkx21 = step*(-ka2*x2(i)+kb2*I(i));
	rkx31 = step*(-ka3*x3(i)+kb3*I(i));
	
	rkx12 = step*(-ka1*(x1(i)+rkx11/2)+kb1*I(i));
	rkx22 = step*(-ka2*(x2(i)+rkx21/2)+kb2*I(i));
	rkx32 = step*(-ka3*(x3(i)+rkx31/2)+kb3*I(i));
	
	rkx13 = step*(-ka1*(x1(i)+rkx12/2)+kb1*I(i));
	rkx23 = step*(-ka2*(x2(i)+rkx22/2)+kb2*I(i));
	rkx33 = step*(-ka3*(x3(i)+rkx32/2)+kb3*I(i));
	
	rkx14 = step*(-ka1*(x1(i)+rkx13)+kb1*I(i));
	rkx24 = step*(-ka2*(x2(i)+rkx23)+kb2*I(i));
	rkx34 = step*(-ka3*(x3(i)+rkx33)+kb3*I(i));
	x1(i+1) = x1(i)+(rkx11+2*rkx12+2*rkx13+rkx14)/6;
	x2(i+1) = x2(i)+(rkx21+2*rkx22+2*rkx23+rkx24)/6;
	x3(i+1) = x3(i)+(rkx31+2*rkx32+2*rkx33+rkx34)/6;
	
    time(i) = t(i)/60; %Unit of time is hour
    gbasul(i) = 5;     %g_based
    g_high(i) = 8;     
    g_low(i) = 4;     
    
    u(i+1) = u(i);
end

figure(1)
hold on;
title("Eat Quantity");
xlabel('Time (h)');
ylabel('D (g CHO/min)');
plot(time,d_cho(1:length(t)));
hold off;
legend('Meal'); 

figure(2)
hold on;
title("Insulin administration");
xlabel('Time (h)'); 
ylabel('Insulin infusion rate (mU/min)');
plot(time,u(1:length(t)));
hold off;
legend('Insulin infusion rate');

figure(3)
hold on;
title("Plasma glucose concentration");
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
