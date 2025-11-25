clear all; close all; clc;

t = [0;
1800.495546;
3600.991092;
5401.486637;
7201.982183;
9002.477729;
10802.97327;
12603.46882;
14403.96437;
16204.45991;
18004.95546;
19805.451;
21605.94655;
23406.4421;
25206.93764];

xi1 = (1.0e-08)*[0.712248399049116;
  -0.000136546650917;
   0.000000009665067;
  -0.000000000000221];

pos0 = 309.291*10^-6;

dt =  0.01;
t = 0:dt:floor(t(end));

SV0 = [7*10^-6];

% Numerically integrate the displacement
[ts,SV] = ode45(@(t,SV) SingleCellDyn(t,SV,xi1),t,SV0);
figure(20)
hold on
plot(ts,SV(:,1),'b') 
legend('Not Drugged x_{d} (m)','Drugged x_{d} (m)','Interpreter','latex')


function rhs = SingleCellDyn(ts,SV,xi1)
% a = [1,ts,ts.*ts,ts.*ts.*ts,SV(1),SV(1).*SV(1),SV(1).*SV(1).*SV(1),SV(1).*SV(1).*SV(1).*SV(1)];
a = [1,ts,ts.*ts,ts.*ts.*ts];

rhs = (a*xi1)';
end
