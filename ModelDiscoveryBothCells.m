% all

%% |==============================================================================================|%
  %|      Filename: ModelDiscoverySingleCellDrugged.m                                             |%
  %|                                                                                              |%
  %|      Author  : Eric Havenhill, PhD Student                                                   |%
  %|                Cellular Engineering and Mechanobiology Lab                                   |%
  %|                Colorado State University                                                     |%
  %|                                                                                              |%
  %|      Purpose : This code performs a DMD on one cell's trajectory to get an ODE               |%
  %|                                                                                              |%
  %|      Notes   : * Trajectory data obtained from TrackMate                                     |%
%% |==============================================================================================|%

clear all; close all; clc;

format('long')

x = (10^-6)*[271.17;
    272.21;
    274.28;
    276.35;
    283.59;
    288.77;
    291.87;
    292.91;
    294.98;
    297.05;
    294.98;
    294.98;
    294.98];

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

% sol in the workspace
%u_94 = sol(:,100-6);

t = linspace(0,t(end),length(x));
dt = t(2);

%% Model Discovery
% State variables

figure(),plot(t,(10^6)*x,'Linewidth',[2])
xlabel('Time (seconds)'),ylabel('Position $(\mu m)$ : x(t) from TrackMate','Interpreter', 'latex')
title('Position vs. Time')
% axis([0 25206.93764 -5 25])

% Compute the derivatives from this data
% This is the b in Ax=b
n = length(t);
for j = 2:n-1
    xdot(j-1) = (x(j+1)-x(j-1))/(2*dt);
end

% Note: We are only calculating the derivatives on the internal
%       points and will discard the external points. 
xs = x(2:n-1);
ts = t(2:n-1);

%A = [ones(length(xs),1),ts',(ts').*(ts'),(ts').*(ts').*(ts'),xs,xs.*xs,xs.*xs.*xs,xs.*xs.*xs.*xs];
%A = [ones(length(xs),1),xs,(xs).*(xs),(xs).*(xs).*(xs)];
%A = [ones(length(xs),1)];
%A = [ones(length(xs),1),ts',(ts').*(ts')];

%A = [ones(length(xs),1)];
%A = [ones(length(xs),1),ts',(ts').*(ts')];
A = [ones(length(xs),1),ts',(ts').*(ts'),(ts').*(ts').*(ts')];

%% Over-determined
% This system of equations can be easily solved by the \
xi1 = A\xdot.';

figure(2)
bar(xi1),xlabel('Terms in the A-matrix'),ylabel('$\dot{x}(t)$','Interpreter', 'latex')
% xticklabels({'1','t','t^{2}','t^{3}','u^{94}','(u_{94})^{2}','(u_{94})^{3}','(u_{94})^{4}'}) %From the A-matrix
xticklabels({'1','t','t^2','t^{3}'}) %From the A-matrix

% This system of equations can be easily solved by the pseudo-inverse
% xi1 = pinv(A)*u_6dot.';

% Uncomment the following in order to use the Lasso algorithm to promote sparsity
% xi1 = lasso(A,u_6dot.','Lambda',0.01);
%% ODE Solver for the SingleCellDyn
pos0 = 271.17*10^-6;

dt =  0.01;
ts = 0:dt:floor(t(end));

SV0 = [x(1)];

% Numerically integrate the displacement
[ts,SV] = ode45(@(ts,SV) SingleCellDyn(ts,SV,xi1),ts,SV0);
figure(10)
plot(ts,(10^6)*SV(:,1),'Linewidth',[2]') 
xlabel('Time (seconds)'),ylabel('Position $(\mu m)$ : x(t) from DMD','Interpreter', 'latex')
title('Position vs. Time')

% figure(4)
% vel = xi1(1);
% %vel = xi1(1) + xi1(2)*ts + xi1(3)*ts.*ts;
% %vel = xi1(1) + xi1(2)*ts + xi1(3)*ts.*ts + xi1(4)*ts.^3;
% %vel = [xi1']*[ones(length(ts));ts;ts.*ts;ts.^3];
% plot(ts,vel)
% xlabel('Time (seconds)'),ylabel('Velocity (m/s): v_{94}(t)')
% title('Velocity vs. Time')

t_control = [0;
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

x_control = 10^(-6)*[331.0615;
  340.4371;
  348.3310;
  352.2970;
  357.4406;
  359.6030;
  358.9961;
  363.1478;
  362.8268;
  362.2703;
  367.1433;
  369.2125;
  370.0408;
  369.9842;
  369.3377];

% sol in the workspace
%u_94 = sol(:,100-6);

t_control = linspace(0,t_control(end),length(x_control));
dt = t_control(2);

%% Model Discovery
% State variables

figure(),plot(t_control,(10^6)*x_control,'Linewidth',[2])
xlabel('Time (seconds)'),ylabel('Position $(\mu m)$ : x(t) from TrackMate','Interpreter', 'latex')
title('Position vs. Time')
% axis([0 25206.93764 -5 25])

% Compute the derivatives from this data
% This is the b in Ax=b
n = length(t_control);
for j = 2:n-1
    xdot_control(j-1) = (x_control(j+1)-x_control(j-1))/(2*dt);
end

% Note: We are only calculating the derivatives on the internal
%       points and will discard the external points. 
xs_control = x_control(2:n-1);
ts_control = t_control(2:n-1);

%A = [ones(length(xs),1),ts',(ts').*(ts'),(ts').*(ts').*(ts'),xs,xs.*xs,xs.*xs.*xs,xs.*xs.*xs.*xs];
%A = [ones(length(xs),1),xs,(xs).*(xs),(xs).*(xs).*(xs)];
%A = [ones(length(xs),1)];
%A = [ones(length(xs),1),ts',(ts').*(ts')];

%A = [ones(length(xs),1)];
%A = [ones(length(xs_control),1),ts_control',(ts_control').*(ts_control')];
A = [ones(length(xs_control),1),ts_control',(ts_control').*(ts_control'),(ts_control').*(ts_control').*(ts_control')];

%% Over-determined
% This system of equations can be easily solved by the \
xi1_control = A\xdot_control.';

figure(2)
bar(xi1),xlabel('Terms in the A-matrix'),ylabel('$\dot{x}(t)$','Interpreter', 'latex')
% xticklabels({'1','t','t^{2}','t^{3}','u^{94}','(u_{94})^{2}','(u_{94})^{3}','(u_{94})^{4}'}) %From the A-matrix
xticklabels({'1','t','t^2','t^{3}'}) %From the A-matrix

% This system of equations can be easily solved by the pseudo-inverse
% xi1 = pinv(A)*u_6dot.';

% Uncomment the following in order to use the Lasso algorithm to promote sparsity
% xi1 = lasso(A,u_6dot.','Lambda',0.01);
%% ODE Solver for the SingleCellDyn
pos0 = 271.17*10^-6;

dt =  0.01;
ts_control = 0:dt:floor(t_control(end));

SV0_control = [x_control(1)];

% Set figure dimensions and font size
fig_width = 600;
fig_height = 500;
font_size = 10; % Set consistent font size

% Numerically integrate the displacement
[ts_control,SV_control] = ode45(@(ts_control,SV) SingleCellDynControl(ts_control,SV0_control,xi1_control),ts_control,SV0_control);
figure(10)
plot(ts_control,(10^6)*SV_control(:,1)-(10^6)*SV_control(1,1)+7,'Linewidth',[2]') 
xlabel('Time (seconds)'),ylabel('Position','Interpreter', 'latex')
title('Position vs. Time')
hold on 

plot(ts,(10^6)*SV(:,1)-(10^6)*SV(1,1)+7,'Linewidth',[2]') 
xlabel('Time (seconds)'),ylabel('Position','Interpreter', 'latex')
title('Position vs. Time')
hold on;

% Create legend with proper sizing
leg = legend('Not Drugged $x_{d}$ ($\mu$m)','Drugged $x_{d}$ ($\mu$m)','Interpreter','latex','FontSize', font_size, 'Location', 'northwest');

% Adjust legend box size and position
leg.Position(3) = leg.Position(3) * 2.05; % Increase width by 30%
leg.Position(4) = leg.Position(4) * 1; % Increase height by 20%
leg.Position(1) = 0.149; % Set x-position (0 = left edge, 1 = right edge)
leg.Position(2) = 0.83; % Set y-position (0 = bottom edge, 1 = top edge)

% Set axes font size and figure dimensions
ax = gca;
ax.FontSize = 15;
set(gcf, 'Position', [100, 100, fig_width, fig_height]);

% 
% figure(4)
% vel = xi1(1);
% %vel = xi1(1) + xi1(2)*ts + xi1(3)*ts.*ts;
% %vel = xi1(1) + xi1(2)*ts + xi1(3)*ts.*ts + xi1(4)*ts.^3;
% %vel = [xi1']*[ones(length(ts));ts;ts.*ts;ts.^3];
% plot(ts,vel)
% xlabel('Time (seconds)'),ylabel('Velocity (m/s): v_{94}(t)')
% title('Velocity vs. Time')

function rhs = SingleCellDyn(ts,SV,xi1)
%a = [1,ts,ts.*ts,ts.*ts.*ts,SV(1),SV(1).*SV(1),SV(1).*SV(1).*SV(1),SV(1).*SV(1).*SV(1).*SV(1)];

%a = [1];
%a = [1,ts,ts.^2];
a = [1,ts,ts.^2,ts.^3];
%a = [1];

rhs = (a*xi1)';
end

function rhs = SingleCellDynControl(ts,SV,xi1)
%a = [1,ts,ts.*ts,ts.*ts.*ts,SV(1),SV(1).*SV(1),SV(1).*SV(1).*SV(1),SV(1).*SV(1).*SV(1).*SV(1)];

%a = [1];
%a = [1,ts,ts.^2];
a = [1,ts,ts.^2,ts.^3];
%a = [1];

rhs = (a*xi1)';
end