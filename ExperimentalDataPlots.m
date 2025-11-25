clear all; close all; clc;

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

x_control = [  331.0615;
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

plot(t_control,x_control-x_control(1,1)+7,'b')
xlabel('Time (seconds)'),ylabel('Position ($\mu$ m)','Interpreter','Latex')
title('Position vs. Time')
hold on

x_drugged = [271.17;
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

length(x_drugged)

t_drugged = [0;
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
    23406.442];

length(t_drugged)

    % Set figure dimensions (width x height in pixels)
    fig_width = 600;
    fig_height = 500;
    set(gcf, 'Position', [100, 100, fig_width, fig_height]);
    % Set font size (customizable)
    font_size = 10; % Change this value as needed

plot(t_drugged,x_drugged-x_drugged(1,1)+7,'r')
xlabel('Time (seconds)'),ylabel('Position ($\mu$ m)','Interpreter','Latex')
title('Position vs. Time')
leg = legend('Not Drugged $x_{d}$ (m)','Drugged $x_{d}$ (m)','Interpreter','latex','FontSize', font_size, 'Location', 'northwest')
    % Adjust legend box size and position - now this will work
    leg.Position(3) = leg.Position(3) * 1.2; % Increase width by 300%
    leg.Position(1) = leg.Position(1) + 0; % Move right by 0.05 (adjust as needed)
    leg.Position(2) = leg.Position(1) + 0.675; % Move right by 0.05 (adjust as needed)
    % Set axes font size and y-axis limits
    ax = gca;
    ax.FontSize = font_size;
    ax = gca; ax.FontSize = 15;

