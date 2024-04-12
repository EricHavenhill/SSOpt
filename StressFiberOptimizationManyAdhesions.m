%% |==============================================================================================|%
  %|      Filename: StressFiberOptimizationManyAdhesions.m                                        |%
  %|                                                                                              |%
  %|      Author  : Eric Havenhill, PhD Student                                                   |%
  %|                Cellular Engineering and Mechanobiology Lab                                   |%
  %|                Colorado State University                                                     |%
  %|                                                                                              |%
  %|      Purpose : This code performs a structural optimization (sizing) on a cell with 7        |%
  %                 elements serving as stress fibers                                             |%
  %|                                                                                              |%
  %|      Notes   : * The EOM are obtained from the FEM                                           |%
  %|                * Static dynamics are used                                                    |%
  %|                * Initial position (x,y,z) -> assumed                                         |%
  %|                * Unit System  = kg,m,s                                                       |%
%% |==============================================================================================|%

clc; clear all; close all;
    
global E Fa sigma_max L_mid theta L Ne Nn

Ne = 7;
Nn = Ne +1;

% Small Elasticity
% E = 0.5*10^6;

% Mid. Elasticity 
% E = 1.45*10^6;
% 
% % Large Elasticity
E = 2.4*10^6;

% Fa = 1*10^-6;
Fa = 120*10^-9;
sigma_max=12*10^6;

% For the non-protruding (symmetric) case:
% theta = (pi/180)*linspace(15,180-15,Ne)';

% For the protruding (offset) case:
theta = (pi/180)*linspace(90,180-15,Ne)';

L_mid = 26*10^-6;

for n = 1 : Ne
    L(1,n) = L_mid/sin(theta(n,1));
end

%% Initial guess for the values
% X = [A1,...,A11,ux,uy]
X0 = [(10^-6)*ones(Ne,1);(10^-6)*ones(2,1)]; 

%% Use fmincon
X = fmincon(@MyNLObj,X0,[],[],[],[],[],[],...
    @MyNLCon)

for i = 1: Ne
    FiberForce(i) = X(i)*E*sqrt(X(end-1)^2 + X(end)^2)/L(i)
end

%% Write the cost function
function J = MyNLObj(X)
    global E Fa sigma_max L_mid theta L Ne Nn
    % Enter the cost function, below. Note: A1 = X(1,1) & A2 = X(2,1) & etc. 
    J = L*X(1:Ne,1);
end

%% Write the nonlinear constraint function
function [g,geq] = MyNLCon(X)
    global E Fa sigma_max L_mid theta L Ne Nn
    
    % Enter the constraint function, below
    % The below function is the inequality 
    % constraint
    
    for i=1:Ne
        % Stress constraint
        g1(i,1) = -sigma_max - (E/L(i))*...
            (cos(theta(i))*X(end-1,1) + ...
            sin(theta(i))*X(end,1));
        g2(i,1) = -sigma_max + (E/L(i))*...
            (cos(theta(i))*X(end-1,1) ...
        + sin(theta(i))*X(end,1));
        % Buckling constraint
        g3(i,1) = cos(theta(i))*X(end-1,1) ...
            + sin(theta(i))*X(end,1) - ...
            (pi*X(i,1)/(4*L(i)));
    end
    
    % Constrain the areas to be bounded:
    % 0.75*10^(-6) <= A_p <= 4.5*10^(-6)
    % Note that weights are placed on both 
    % of these terms to more strongly enforce
    % these constraints
    g4 = 50*(0.75*10^(-6)-X(1:end-2,1));
    g5 = 50*(-4.5*10^(-6)+X(1:end-2,1));
    
    g  = [g1;g2;g3;g4;g5];
    
    % Formulate the equality constraint 
    % based on the FEM formulation of the
    % equilibrium equaitons
    
    k11 = 0;
    k12 = 0;
    for i = 1:Ne
        k11 = k11 + X(i,1)*E*...
            cos(theta(i))^(2)/L(i);
        k12 = k12 + X(i,1)*E*...
            cos(theta(i))^(2)*...
            sin(theta(i))^(2)/L(i);
    end
    
    k21=k12;
    k22=-k11; 
    
    geq1 = Fa - k11*X(end-1,1) - k12*X(end,1);
    geq2 = 0 - k21*X(end-1,1) + k22*X(end,1);
    
    % Constrain the y-disp. to be = 0
    geq3 = X(end,1);
    
    geq = [geq1;geq2;geq3];
end
