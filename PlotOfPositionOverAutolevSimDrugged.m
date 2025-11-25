clear all; close all; clc;

t = [0;
    1800.496;
    3600.991;
    5401.487;
    7201.982;
    9002.478;
    10802.973;
    12603.469;
    14403.964;
    16204.46;
    18004.955;
    19805.451;
    21605.947;
    23406.442];

xi1 =(1E-8)*[-0.139666289694561;
            0.000119586902340;
            -0.000000011216122;
            0.000000000000276];

dt =  0.01;
t = 0:dt:floor(t(end));

SV0 = [7*10^-6];

% Numerically integrate the displacement
[ts,SV] = ode45(@(t,SV) SingleCellDyn(t,SV,xi1),t,SV0);
hold on
figure(20)
hold on
plot(ts,SV(:,1),'r') 
legend('Not Drugged $x_{d}$ (m)','Drugged $x_{d}$ (m)','Interpreter','latex')
title('Position vs. Time')
xlabel('Time (seconds)')
ylabel('Position (m)','Interpreter','latex')
box on


function rhs = SingleCellDyn(ts,SV,xi1)
% a = [1,ts,ts.*ts,ts.*ts.*ts,SV(1),SV(1).*SV(1),SV(1).*SV(1).*SV(1),SV(1).*SV(1).*SV(1).*SV(1)];
a = [1,ts,ts.*ts,ts.*ts.*ts];

rhs = (a*xi1)';
end