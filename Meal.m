clc
clear
close all
a=0;
b=1500;
step = 0.01;
t=a:step:b; 
%parameter for Meals
tau_d = 40;
Ag = 0.8;
Dm1(1) = 0;
Dm2(1) = 0;
%Molecular weight of glouse
Mwg = 180;
for j = 1 : length(t)
    %equal  d_cho(j)=10*heaviside(t(j)-1)-10*heaviside(t(j)-1.2);
    if t(j)>=60 && t(j)<=70
        d_cho(j)=10;
    else
		d_cho(j)=0;
    end
    D_meal(j) = 1000*d_cho(j)/Mwg;
    
	k11 = Ag*D_meal(j)-Dm1(j)/tau_d;
    k21 = Dm1(j)/tau_d-Dm2(j)/tau_d;
    
    k12 = Ag*D_meal(j)-(Dm1(j)+step/2*k11)/tau_d;
    k22 = (Dm1(j)+step/2*k11)/tau_d-(Dm2(j)+step/2*k21)/tau_d;
    
    k13= Ag*D_meal(j)-(Dm1(j)+step/2*k12)/tau_d;
    k23= (Dm1(j)+step/2*k12)/tau_d-(Dm2(j)+step/2*k22)/tau_d;
    
    k14= Ag*D_meal(j)-(Dm1(j)+step*k13)/tau_d;
    k24= (Dm1(j)+step*k13)/tau_d-(Dm2(j)+step*k23)/tau_d;
    
	Dm1(j+1) = Dm1(j)+step*(k11+2*k12+2*k13+k14)/6;
	Dm2(j+1) = Dm2(j)+step*(k21+2*k22+2*k23+k24)/6;
    
	Ug(j) = Dm2(j)/tau_d;
end

figure(1)
hold on;
xlabel('Time (h)'); 
ylabel('D (g CHO/min)');
plot(t,d_cho(1:length(t)));
hold off;
legend('Meal', 'Meal2'); 

figure(2)
hold on;
xlabel('Time (h)'); 
ylabel('Ug (g CHO/min)');
plot(t,Ug(1:length(t)));
hold off;
legend('Meal', 'Meal2'); 

