clc
clear 
close All
%% Sec.1 Initial Values

r = 8; %Compression Ratio
gamma= 1.4; %ratio of specific heats
ro=1.225; %air density 
N=5000; % RpM
V_d=0.001998/4; %Displacement (cc)
V_bdc=(r/(r-1))*V_d;
V_tdc=V_bdc/r;
s = 0.086; %Stroke (m) 
a1=s/2; % crankshaft radius
l_cr= 0.153; %connecting rod length(m)
theta=-360:1:720; %crank angle position
epsilon = a1/l_cr;
b=((4*V_d)/(s*pi))^(0.5); %bore;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sec.2 calculating volume & area , and figures

V=((pi/4)*(((b.^2)*s)/(r-1)))+((pi/4)*(b.^2))*(l_cr+a1-a1*cosd(theta)-((l_cr.^2)-(a1*sind(theta)).^2).^0.5);
A=((pi/2)*(b.^2))+((pi*b*s*0.5)*(l_cr+a1-a1*cosd(theta)-((l_cr.^2)-(a1*sind(theta)).^2).^0.5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(1)
plot(theta,1000000*V,'linewidth',2,'color',[1 0 0]);
hold on;
plot(theta,1000000*V_bdc*ones(size(theta)),'-.','linewidth',1.5,'color',[0 0 1]);
hold on;
plot(theta,1000000*V_tdc*ones(size(theta)),'-.','linewidth',1.5,'color',[0 0 1]);

figure(2)
plot(theta,A,'linewidth',3,'color',[1 0.46 0.43]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sec.3 Velocity and its figure
velocity=((pi/2*sind(theta)).*((1+(epsilon*cosd(theta))./(1 - epsilon^2*sind(theta).^2).^(1/2))))*2*pi*N*s/60; % m/s
figure(3)
plot(theta,velocity,'linewidth',3,'color',[0 1 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sec.4 Calculating Burn Fraction

a = 4.5;%2-7
n = 3;%1-5
theta_s1=-16;
diuration_theta = 60;  %40-80 degree
f=(-1.*(theta<theta_s1))+(1.-exp(-a.*((theta-theta_s1)./diuration_theta).^n).*(theta>=theta_s1 ));
figure(4)
plot(theta,f,'r','linewidth',2);
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sec.5 Calculating Unburned/Burned Area,Mass,Temperature,Volume & Motoring Pressure,Pressure, and Work
%preparing Variables
 D_qw(1:360)=zeros;
 D_qw2(1:360)=zeros;
 M_F(1:360)=zeros;
 DQ2(1:360)=zeros;
 m_u(1:360)=zeros;
 m_b(1:360)=zeros;
 V_u(1:360)=zeros; 
 flux(1:360)=zeros; 
 rho(1:360)=zeros; 
 V(1:360)=zeros;
 DT(1:360)=zeros;
 Pm(1:360)=zeros;
 DP(1:360)=zeros;
 DT_u(1:360)=zeros;
 P(1:360)=zeros;
 T(1:360)=zeros;
 W(1:360)=zeros;
 V_b(1:360)=zeros;
 T_b(1:360)=zeros;
 T_u(1:360)=zeros;
 DQ_w(1:360)=zeros;
 DQ(1:360)=zeros;
 X(1:360)=zeros;
 DX(1:360)=zeros;
 h_g(1:360)=zeros;
 Q(1:360)=zeros;
 Q2(1:360)=zeros;
 D_W(1:360)=zeros;
 DV(1:360) = zeros;
 A(1:360)=zeros;
 A_u(1:360)=zeros;
 U(1:360)=zeros;
 A_b(1:360)=zeros;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_a = ro*V_d; %Air Mass
m_f = m_a/(1.3*15.09); %Fuel Mass
m_c = m_a+m_f; %Mass
Ne=2000/60;
Up=2*Ne*s; % mean piston speed
P_atm = 101000; %atm Pressure
P_BDC = P_atm;
AF_ratio_ac=1.3*14.6;
T_BDC=300; %Temperature
T(1) = T_BDC; 
IVC = -131; %intake valve close
EVO = 125; %exhaust valve open
lambda=1;
N_cyl=4;
theta_b = 60 ; 
theta_0 = -16; %Start of combustion
Tw=360;
theta_f = theta_0+theta_b; % End of combustion
LHV = 44.6e6;
R_air = 287;
R=R_air/1000;
eta_combmax=0.95;
eta_comb=eta_combmax*(-1.6082+4.6509*lambda-2.0764*lambda^2); %combustion efficiency 
V(1)=V_bdc;
X(1) = 0;
DX(1) = 0;
A(1)=0.0126;
DV(1) = 0;
L = 0.086; 
B = 0.086; 
A_p = (pi/4)*b^2;
T_BDC=300;
T(1) = T_BDC;
P(1)=P_BDC;
V_u(1) = V(1);
V_b(1)=0;
A_ch = A_p;
T_u(1)=T_BDC;
m_u(1)=m_c;
m_b(1)=0;
theta=(-180:179);
rho(1) = P(1)/(R_air*T(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:360
rho(i) = P(i-1)/(R_air*T(i-1));
V(i)=((pi/4)*(((b.^2)*s)/(r-1)))+((pi/4)*(b.^2))*(l_cr+a1-a1*cosd(theta(i))-((l_cr.^2)-(a1*sind(theta(i))).^2).^0.5);
DV(i) = V(i)-V(i-1);
A(i)=((pi/2)*(b.^2))+((pi*b*s*0.5)*(l_cr+a1-a1*cosd(theta(i))-((l_cr.^2)-(a1*sind(theta(i))).^2).^0.5));
X(i)=(-1.*(theta(i)<theta_0))+(1.-exp(-a.*((theta(i)-theta_0)./theta_b).^n).*(theta(i)>=theta_0 ));
DX(i) = X(i) - X(i-1);
M_F(i) =(( V(164)*rho(164)/(lambda*AF_ratio_ac)).*(theta(i) < theta_f));
A_u(i)=A(i)*(1-sqrt(X(i)));
A_b(i)=A(i)*(X(i)/sqrt(X(i)));

D_qw(i)=((h_g(i)*A(i)/(360*Ne))*(T(i-1)-Tw));
D_qw2(i)=(((h_g(i)*A_u(i)/Ne)/4)*(T_u(i-1)-Tw))+(((h_g(i)*A_b(i)/6*Ne)/4)*(T_b(i-1)-Tw));
DQ(i) = eta_comb*LHV*M_F(i)*DX(i)-D_qw(i);
DQ2(i) = eta_comb*LHV*M_F(i)*DX(i)-D_qw2(i);
Q(i) = Q(i-1)+DQ(i);
Q_loss= (N_cyl*Q(i)*Ne/2)/1000;

Q2(i) = Q2(i-1)+DQ2(i);
if IVC< theta(i)
 DT(i)=T(i-1)*(gamma-1)*((1/(P(i-1)*V(i-1)))*DQ(i)-(1/V(i-1))*DV(i));
 DP(i)=(((-P(i-1)/V(i-1))*DV(i))+((P(i-1)/T(i-1))*DT(i)));


end
  P(i) = P(i-1)+DP(i);

 if EVO < theta(i)
 if P(i)<=P_atm
 P(i)=P_atm;
 end
 end
 

 if theta(i)< IVC
     P(i) = P_atm;
 end
  T(i) = T(i-1)+DT(i);
  
  Pm(i)=P(50)*((V(50)./V(i))^(gamma));
 if theta(i)<theta(50)
     Pm(i)=P_BDC;
 end
 m_b(i) = m_b(i-1)+DX(i)*m_c; 
 m_u(i) = m_u(i-1)-DX(i)*m_c;
 if theta(i)<=theta_0
 %V_u(i)=N_cyl*V(i);
  V_u(i)=V(i);

 end
 if theta(i)>theta_0
 V_u(i)=(((m_u(i)*V_u(i-1)))/m_u(i-1))*(P(i)/P(i-1))^(-1/gamma);
 end
 %V_b(i)=(N_cyl*V(i))-V_u(i);
  V_b(i)=(V(i))-V_u(i);

 if V_b(i)<0
 V_b(i)=0;
 end
  T_u(i)=P(i)*V_u(i)/(m_u(i)*R*1000);
 if theta(i) <= theta_0
 T_b(i)=0;
 end
 if theta(i)>theta_0
 T_b(i)=P(i)*V_b(i)/(m_b(i)*R*1000);
 end

 DT_u(i)=T_u(i)-T_u(i-1); 
 U(i)=(((2.28*Up)+(0.00324*(T(50))*(V_d/V(50))*((P(i)-Pm(i))/P(50)))).*(theta(i)>=IVC & theta(i)<=EVO))+((6.18*Up).*(theta(i)<IVC & theta(i)>EVO));
 h_g(i)=3.26*((P(i))^(0.8))*((U(i))^0.8)*(b^(-0.2))*((T(i))^(-0.55));
 W(i) = W(i-1)+(P(i)-P_atm)*DV(i);
 D_W(i)=(N_cyl*W(i)*Ne/2)/1000;
 flux(i) = h_g(i) * (T(i) - Tw) ;
   
  %if theta(i)> theta_f+15
  %   T(i)= T_b(i);
  %end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sec.6 Some Values for Plots

y_A=linspace(0,0.02,100);
y_F=linspace(0,1.5,100);
y_M=linspace(0,7e-04,100);
y_T=linspace(-1000,3000,1000);
y_V=linspace(0,700,1000);
y_P=linspace(0,45,1000); 

%% Sec.7 Figures
figure(1)
title('\it Volume(\bf\theta\circ)','fontsize',18,'color',[1 1 0]);
text(-260,100,'V_t_d_c=71.35');
text(-60,600,'V_b_d_c=570.8');
set(gca,'color',[0.5 0.5 0.5],'Xlim',[-360 360],'Xcolor',[1 1 0],'Ylim',[0 700],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it Crank Angle \bf\theta\circ','fontsize', 11, 'color',[1 1 0]);
ylabel('\it Each Cylinder Volume [cc]','fontsize', 11, 'color',[1 1 0]);
legend( 'Volume');
grid on;

%%%%%%%%%

figure(2)
title('\it Area(\bf\theta\circ)','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.5 0.5 0.5],'Xlim',[-360 360],'Xcolor',[1 1 0],'Ylim',[112e-04 130e-04],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it Crank Angle \bf\theta\circ','fontsize', 11, 'color',[1 1 0]);
ylabel('\it Cylinder Area [m^2]','fontsize', 11, 'color',[1 1 0]);
legend('Area');
grid on;

%%%%%%%%%
figure(3)
title('\it Velocity(\bf\theta\circ)','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.5 0.5 0.5],'Xlim',[-360 360],'Xcolor',[1 1 0],'Ylim',[-80 80],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it Crank Angle \bf\theta\circ','fontsize', 11, 'color',[1 1 0]);
ylabel('\it Piston Speed [m/s]','fontsize', 11, 'color',[1 1 0]);
legend('Speed');
grid on;

figure(4)
hold on;
plot(theta_0*ones(size(y_F)),y_F,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
hold on;
plot(theta_f*ones(size(y_F)),y_F,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
title('\it Burn Fraction(\bf\theta\circ)','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.5 0.5 0.5],'Xlim',[-120 120],'Xcolor',[1 1 0],'Ylim',[0 1.5],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it Crank Angle \bf\theta\circ','fontsize', 11, 'color',[1 1 0]);
ylabel('\it Cumulative Burn Fraction','fontsize', 11, 'color',[1 1 0]);
legend('Cylinder 1 Burn Fraction' );
grid on;

%%%%%%%%%

figure(5)
plot(theta,A_u,'linewidth',2,'color',[0.6 0 0.6]);
hold on;
plot(theta,A_b,'linewidth',2,'color',[0 0.7 0.7]);
hold on;
plot(theta_0*ones(size(y_A)),y_A,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
hold on;
plot(theta_f*ones(size(y_A)),y_A,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
title('\it Unburned & Burned Area(\bf\theta\circ)','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.5 0.5 0.5],'Xlim',[-178 180],'Xcolor',[1 1 0],'Ylim',[0 0.014],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it Crank Angle \bf\theta\circ','fontsize', 11, 'color',[1 1 0]);
ylabel('\it Unburned Area [m^2]','fontsize', 11, 'color',[1 1 0]);
legend('A\it_u_n_b_u_r_n_e_d','A\it_b_u_r_n_e_d');
grid on;

%%%%%%%%%

figure(6)
plot(theta,m_u*1000,'linewidth',3,'color',[1 0.5 0.5]);
hold on;
plot(theta,m_b*1000,'linewidth',3,'color',[0.5 1 0.5]);
hold on;
plot(theta_0*ones(size(y_A)),y_A*1000,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
hold on;
plot(theta_f*ones(size(y_A)),y_A*1000,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
title('\it Unburned & Burned Mass(\bf\theta\circ)','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.5 0.5 0.5],'Xlim',[-180 180],'Xcolor',[1 1 0],'Ylim',[0 2],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it Crank Angle \bf\theta\circ','fontsize', 11, 'color',[1 1 0]);
ylabel('\it Mass [gr]','fontsize', 11, 'color',[1 1 0]);
legend('m\it_u_n_b_u_r_n_e_d','m\it_b_u_r_n_e_d');
grid on;

%%%%%%%%%

figure(7)
plot(theta(2:226),T_u(2:226),'linewidth',3,'color',[0.5 1 0.5]);
hold on;
plot(theta(165:360),T_b(165:360),'linewidth',3,'color',[0.5 0.5 1]);
hold on;
plot(theta,T,'-.','linewidth',3,'color',[1 0 0]);
hold on;
plot(theta_0*ones(size(y_T)),y_T*1000,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
hold on;
plot(theta_f*ones(size(y_T)),y_T*1000,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
title('\it Temperature(\bf\theta\circ)','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.5 0.5 0.5],'Xlim',[-180 180],'Xcolor',[1 1 0],'Ylim',[0 3000],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it Crank Angle \bf\theta\circ','fontsize', 11, 'color',[1 1 0]);
ylabel('\it Temperature \bf K \circ','fontsize', 11, 'color',[1 1 0]);
legend('T\it_u_n_b_u_r_n_e_d' , 'T\it_b_u_r_n_e_d' , 'T' );
%legend('T\it_u_n_b_u_r_n_e_d' , 'T\it_b_u_r_n_e_d' );

grid on;

%%%%%%%%%

figure(8)
plot(theta,1000000*V_u,'linewidth',3,'color',[0.25 1 0.25]);
hold on;
plot(theta,1000000*V_b,'linewidth',3,'color',[0.25 0.25 1]);
hold on;
plot(theta_0*ones(size(y_V)),y_V*1000,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
hold on;
plot(theta_f*ones(size(y_V)),y_V*1000,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
title('\it Unburned & Burned Volume(\bf\theta\circ)','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.5 0.5 0.5],'Xlim',[-180 180],'Xcolor',[1 1 0],'Ylim',[0 700],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it Crank Angle \bf\theta\circ','fontsize', 11, 'color',[1 1 0]);
ylabel('\it Volume [cc] ','fontsize', 11, 'color',[1 1 0]);
legend('V\it_u_n_b_u_r_n_e_d' , 'V\it_b_u_r_n_e_d');
grid on;

%%%%%%%%%

figure(9)
plot(theta,P/100000 ,'linewidth',2,'color',[1 0.15 0.15]);
hold on;
plot(theta_0*ones(size(y_P)),y_P,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
hold on;
plot(theta_f*ones(size(y_P)),y_P,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
hold on;
plot(0*ones(size(y_P)),y_P,'-.','linewidth',1,'color',[1 1 1]);
title('\it Pressure(\bf\theta\circ)','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.5 0.5 0.5],'Xlim',[-180 180],'Xcolor',[1 1 0],'Ylim',[0 45],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it Crank Angle \bf\theta\circ','fontsize', 11, 'color',[1 1 0]);
ylabel('\it Pressure [bar] ','fontsize', 11, 'color',[1 1 0]);
legend('P' );
grid on;

%%%%%%%%%

figure(10)
plot(theta,W,'linewidth',2,'color',[0.25 1 0.25]);
hold on
title('\it Work(\bf\theta\circ)','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.5 0.5 0.5],'Xlim',[-180 180],'Xcolor',[1 1 0],'Ylim',[-200 800],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it Crank Angle \bf\theta\circ','fontsize', 11, 'color',[1 1 0]);
ylabel('\it Work [J] ','fontsize', 11, 'color',[1 1 0]);
legend('W');
grid on;

%%%%%%%%%

figure(11)
plot(theta,Pm/100000,'linewidth',2,'color',[0.5 1 0.75]);
hold on;
plot(theta_0*ones(size(y_P)),y_P,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
hold on;
plot(theta_f*ones(size(y_P)),y_P,'-.','linewidth',0.9,'color',[0.85 0.85 0.4]);
title('\it Motoring Pressure(\bf\theta\circ)','fontsize',18,'color',[1 1 0]);
set(gca,'color',[0.5 0.5 0.5],'Xlim',[-180 180],'Xcolor',[1 1 0],'Ylim',[0 45],'Ycolor',[1 1 0],'fontsize',11,'linewidth',2);
set(gcf,'color',[0 0 0])
xlabel('\it Crank Angle \bf\theta\circ','fontsize', 11, 'color',[1 1 0]);
ylabel('\it Pressure [bar] ','fontsize', 11, 'color',[1 1 0]);
legend('P_m');
grid on;