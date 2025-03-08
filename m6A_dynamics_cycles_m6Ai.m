

function m6A_dynamics_cycles
% Parameters

clear all
%close all

global b1 b2 d1 d2 alpha C_Dox1 C_Dox2 Vs_1 Vs_2 T_int C_int T_cyc beta_2 beta_1 B_int

% b1 = 0.1;   % Birth rate in non-resistant tumour cells X_1
% b2 = 0.1;   % Birth rate in resistant tumour cells X_2
% d1 = 0.05;  % Death rate in non-resistant tumour cells
% d2 = 0.05;  % Death rate in resistant tumour cells

alpha = 1*10^(-4);    % m6A methylation coeff by m6A writers
beta_1 = 1*10^(-4);     % m6A de-methylation coeff by m6A erasers
beta_2 = 100*beta_1;

C_Dox1 = [0 0.037 0.07802 0.1509 0.2998 0.5989 1.199]; % SW48
Vs_1 = [1.0 0.5853 0.3616 0.28 0.1923 0.1406 0.05147];
C_Dox2 = [ 0 0.044 0.082 0.1525 0.3024 0.5999 1.2];     % TP53 Vehicle
Vs_2 =  [ 1.0 0.9012 0.807 0.7112 0.6213 0.46441 0.2879];

C_int = [.67 .67 .67]; % interval dosing in micro Molar
B_int = [0 0 1]; % interval rate for beta

% Initial Conditions; Total cell conc = 1 micro Molar
x1_count = 95; % micro Molar
x2_count = 5; % micro Molar
Init = [x1_count, x2_count];

% Time Span
T_int = 48; 
T_cyc = 3*T_int;
Tfinal = 3*T_cyc;
tspan = [0, Tfinal];



% Check death rates and birth rates

%c = 0.1*10^(-6); % Dox at 0.1 micro Molar
c = 0.1; % mM
rs_D1 = growth_rate_X1(c);
rs_D2 = growth_rate_X2(c);
t_0 = 3*24;
t_f = 5*24;
Af_1 =.2553; % A/A+S for non resistant X1 SW 48
Af_2 = .0504; % A/A+S for resistant X2    TP53

d_c1 = Af_1*(exp(rs_D1*(t_f)))/( ((exp(t_f*rs_D1)-exp(t_0*rs_D1))/rs_D1)*(1 - Af_1) ); % Death rate calc
b1 = rs_D1 + d_c1; % Birth rate calc
rs_d0_1 = log(2)/35;     % growth rate for X1, with no Dox
Via_c_test1 = exp(rs_D1)/exp(rs_d0_1); % Viability % for X1 to check with plot

d_c2 = Af_2*(exp(rs_D2*(t_f)))/( ((exp(t_f*rs_D2)-exp(t_0*rs_D2))/rs_D2)*(1 - Af_2) ); % Death rate calc
b2 = rs_D2 + d_c2; % Birth rate calc
rs_d0_2 = log(2)/24;    % growth rate for X2, with no Dox
Via_c_test2 = exp(rs_D2)/exp(rs_d0_2); % Viability % for X2 to check with plot
 
% % Viability curve with all Dox concentrations
 
Dox_c = linspace(0, 1.19,100);
rs_Dox_c1 = growth_rate_X1(Dox_c);
rs_Dox_c2 = growth_rate_X2(Dox_c);

Via_1 = exp(rs_Dox_c1*72)/exp(rs_d0_1*72);% 72 hr experiment
Via_2 = exp(rs_Dox_c2*72)/exp(rs_d0_2*72);

% figure;
% plot(Dox_c, Via_1.*100, Dox_c, Via_2.*100);
% hold on
% plot(C_Dox1,Vs_1*100, 'Marker','square','LineStyle','none','Color','blue');
% plot(C_Dox2,Vs_2*100, 'Marker','square','LineStyle','none','Color','green');
% legend('SW48-Dox Model $X_1$', 'TP53-Dox-vehicle Model $X_2$','SW48-Dox Data','TP53-Dox-vehicle Data','interpreter','latex');
% title('Viability ($\%$) for $X1$ and $X2$','Interpreter','latex');
% xlabel('Dox concentration $\mu$ M','Interpreter','latex');
% ylabel('Viability($\%$) ','Interpreter','latex');


figure;
D1 = b1-rs_Dox_c1;D2 = b2-rs_Dox_c2;
D1(D1<0) = 0;
D2(D2<0) = 0;

% plot(Dox_c, D1, Dox_c, D2);
% legend('Unmethylated $X_1$', 'Methylated $X_2$','interpreter','latex');
% title('Death rate for $X1$ and $X2$','Interpreter','latex');
% xlabel('Dox concentration $mM$','Interpreter','latex');
% ylabel('Death rate $D(d)$','Interpreter','latex');
% 
% figure;
% plot(Dox_c, rs_Dox_c1, Dox_c, rs_Dox_c2);
% legend('Unmethylated $X_1$', 'Methylated $X_2$','interpreter','latex');
% title('Growth rate for $X1$ and $X2$','Interpreter','latex');
% xlabel('Dox concentration $mM$','Interpreter','latex');
% ylabel('Growth rate $r_s(d)$','Interpreter','latex');


% Solve ODE
[t, X] = ode45(@(t, X) gene_regulation(t, X), tspan, Init);

%Plot Results
figure(4);
hold on
plot(t, X(:,1),'LineWidth',2,'LineStyle','-');
plot(t, X(:,2),'LineWidth',2,'LineStyle','--');
xlabel('Time [hrs]','interpreter','latex');
ylabel('Cell Concentrations $X1$, $X2$ [$\mu$ M]','interpreter','latex');
legend('$X_1$', '$X_2$','interpreter','latex');
str = ['Dox : [' num2str(C_int(1)) ' ' num2str(C_int(2)) ' ' num2str(C_int(3)) '] $\mu M$ ' '$T_{interval} = $' num2str(T_int) ' hrs ' 'No. of cycles = ' num2str(Tfinal/T_cyc) ];
str2 = [ ' m6Ai : [' num2str(B_int(1)) ' ' num2str(B_int(2)) ' ' num2str(B_int(3)) ']'];
title(str,str2,'interpreter','latex');
hold off

figure(5)
hold on
plot(t, X(:,1)+X(:,2),'LineWidth',2,'LineStyle','-');
ylabel('Total Cell Count [$\mu$ M]','interpreter','latex');
hold off
end

function dy = gene_regulation(t, X)
% Initialize the derivatives
dy = zeros(size(X));
global b1 b2 alpha beta_1 C_int T_int T_cyc beta_2 m6_int B_int

c_drug_1 = C_int(1);
c_drug_2 = C_int(2);
c_drug_3 = C_int(3);


tim = mod(t,T_cyc);

if(tim<T_int)                       % First interval
c = c_drug_1;
beta = B_int(1)*beta_2 + (1-B_int(1))*beta_1;

elseif(tim>T_int && tim<2*T_int)    % Second interval
c = c_drug_2;
beta = B_int(2)*beta_2 + (1-B_int(2))*beta_1;

elseif(tim>2*T_int)                 % Third interval
c = c_drug_3;
beta = B_int(3)*beta_2 + (1-B_int(3))*beta_1;

end

D1 = b1 - growth_rate_X1(c);
D2 = b2 - growth_rate_X2(c);
if(D1<0)
    D1=0;
end
if(D2<0)
    D2=0;
end

dy(1) = (b1 - D1)*X(1) - alpha*X(1) + beta*X(2);
dy(2) = (b2 - D2)*X(2) + alpha*X(1) - beta*X(2);

end

function gr_X1=growth_rate_X1(c_drug)
global b1 C_Dox1 Vs_1
% C = [ 0 0.0416 0.0805 0.1501 0.3006 0.5997 1.199];   % micro Molar
% Vs = [ 1.0 0.5556 0.4729 0.4030 0.2551 0.17107 0.09689]; % viability fraction

C = C_Dox1;Vs = Vs_1;
t_doubling=35;  %Hr
rs_0 = log(2)/t_doubling;
rs_D = (log(Vs) + rs_0*72)/72;
%d_D = b1 - rs_D;
gr_X1 = interp1(C,rs_D,c_drug,'linear','extrap');

end

function gr_X2=growth_rate_X2(c_drug)
global b2 C_Dox2 Vs_2
% C = [ 0 0.044 0.082 0.1525 0.3024 0.5999 1.2];   % micro Molar
% Vs = [ 1.0 0.9012 0.807 0.7112 0.6213 0.46441 0.2879]; % viability fraction

C = C_Dox2;Vs = Vs_2;
t_doubling=24;  %Hr
rs_0 = log(2)/t_doubling;
rs_D = (log(Vs) + rs_0*72)/72;
%d_D = b2 -rs_D;
% gr_X2 = interp1(C*10^(-6),d_D,c_drug);
gr_X2 = interp1(C,rs_D,c_drug,'linear','extrap'); %mM

end
