
%% =======================================================================
% File Name: Project1_DMD_MasterFile.m
% Author: Eric Havenhill
% Purpose: * DMD of cell migration
%               * Develops x(t) for arbitrary cell from DMD
%          * Optimal control + x(t) results
%          * Error analysis
%               * Simulated u(t) vs. actual u(t)
%% =======================================================================

clear all; close all; clc;

global xi right_cell_polyfit right_cell_polyval left_cell_polyfit left_cell_polyval right_cell_polyfit_alt right_cell_polyval_alt
global u_alt t_alt
PolyOrderOfApproxLastCell = 10;

% ---------------------------- Input the Data ----------------------------
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

load x_cell.m;
x = sortrows(x_cell.');


% t = [0;
% 1799.4076;
% 3598.815199;
% 5398.222799;
% 7197.630399;
% 8997.037999;
% 10796.4456;
% 12595.8532;
% 14395.2608;
% 16194.6684;
% 17994.076;
% 19793.4836;
% 21592.8912;
% 23392.2988;
% 25191.7064;
% 26991.114;
% 28790.5216;
% 30589.92919];
% 
% load x_cell_C2_R_all.m
% x = sortrows(x_cell_C2_R_all.');

% ========================= END OF INPUTS ================================

% ---------------------------------- DMD ----------------------------------

for i = 1 : size(x,2)
    for j = 1 : size(x,1)
    u(j,i) = x(j,i) - x(j,1);
    end
end

x_start = x(:,1).'; %starting values of the cells


tint = t;

for i = 1 : size(x,2)
    xint(:,i) = x(1:end,1);
end
%% Data from Experiment
dt=t(2)-t(1); % time snapshots of the data

% The cells must be placed in left-to-right numerical order
%   => 1st is smallest values & last is largest value
u = u.';

[m,n]=size(u);

% Locate the position of the right-most & left-most cell
pos_lastcell = size(u,2);
pos_firstcell = 1;

x_start = x(:,1).'; %starting values of the cells

% ========================================================================
%% Take the derivatives of the data

% Now we produce the RHS of Ax=b using the center difference formula
udot=zeros(m,n-2);
for jj=1:m  % walk through rows (space)
for j=2:n-1  % walk through time
   udot(jj,j-1)=( u(jj,j+1)-u(jj,j-1) )/(2*dt);
end
end

for jj=1:m  % walk through rows (space)
    for j=2:n-1  % walk through time
       udot(jj,j-1)=( u(jj,j+1)-u(jj,j) )/(dt);
    end
end

% Only keep the interrior values -- discard the edge points
udot = udot(2:end-1,1:end);
Udot=reshape((udot.'),[],1);

% derv matrices
dx = x(2)-x(1);

% We are using finite difference matrices to make the derivatives easier
% Finite difference differentiation = matrix multiply
D=zeros(m,m); D2=zeros(m,m);
for j=1:m-1
  D(j,j)=1;
  D(j+1,j)=-1;
%
  D2(j,j+1)=1;
  D2(j+1,j)=1;
  D2(j,j)=-2;
end

D=(1/(dx))*D;
D = D(2:end-1,2:end-1);

D2=D2/(dx^2);
D2=D2(2:end-1,2:end-1);

U=reshape( u(2:end-1,2:end-1).',(n-2)*(m-2) ,1 );

for jj=2:n-1
   ux(jj-1,:)=((D*u(2:end-1,jj)).') + [-u(1,jj)/(dx),zeros(1,m-4),0];  % u_x
   uxx(jj-1,:)=((D2*u(2:end-1,jj)).') + [u(1,jj)/(dx^2),zeros(1,m-4),u(end,jj)/(dx^2)];  % u_xx
end

% Extract the displacement of the RHS cell position
right_cell_polyfit = polyfit(tint(1:end),u(:,pos_lastcell),PolyOrderOfApproxLastCell);
right_cell_polyval = polyval(right_cell_polyfit,tint(1:end));

% Extract the displacement of the LHS cell position
left_cell_polyfit = polyfit(tint(1:end),u(:,pos_firstcell),PolyOrderOfApproxLastCell);
left_cell_polyval = polyval(left_cell_polyfit,tint(1:end));

Ux=reshape(ux,[],1);
Uxx=reshape(uxx,[],1);

% ========================================================================
%% Formulate the library (A matrix) ***
% A = [ones(size(Uxx,1),1) U U.^2 U.^3 U.^4 U.^5 Ux Ux.^2 Ux.^3 Ux.^4 Ux.*U Ux.*U.^2 Ux.*U.^3 Ux.*U.^4 Uxx];
A = [Ux Uxx];
% ========================================================================
%% Solve the system of equations

% If using the \, it implements a QR algorithm -> promoting sparsity
% xi = pinv(A)*Udot;
xi=A\Udot;
% xi=lasso(A,Udot,'Lambda',0.5);

% ========================================================================
%% Plot the results ***
figure
bar(xi)
xlabel('Functions \epsilon A-Matrix'); ylabel('\xi');title('Bar Graph - \xi Coefficient Data');
%xticklabels({'1', 'u', 'u.^2', 'u.^3', 'u.^4', 'u.^5', 'Ux', 'Ux.^2', 'Ux.^3', 'Ux.^4', 'Ux.*u', 'Ux.*u.^2', 'Ux.*u.^3', 'Ux.*u.^4', 'Uxx'})
xticklabels({'Ux','Uxx'})
%  xticklabels({'Ux','Ux.^2','Uxx'})

% ----------- Plot the displacement data of last cell vs. time -----------
figure
plot(tint,u(:,pos_firstcell))
hold on
plot(tint,left_cell_polyval)
hold on
plot(tint,u(:,pos_lastcell))
hold on
plot(tint,right_cell_polyval)
% plot(t(2:end-1),right_cell_polyval)
xlabel('Time (seconds)'); ylabel('u(t): Displacement of Last Cell (\mu m)');
legend('Actual Displacement Data from Trackmate - First Cell','Polynomial Fit - First Cell','Actual Displacement Data from Trackmate - Last Cell','Polynomial Fit - Last Cell');
title('Displacement Data from Trackmate of the Last Cell')

% ----------- Plot the derivative data vs. time -----------
figure
hold on

for i = 1 : width(ux)
    plot(tint(2:end-1),ux') 
    xlabel('Time (seconds)'); ylabel('dudx (\mu m/\mu m)');title('Derivative Data - Interrior Points');
    legendInfo{i} = ['dudx|_x_{',num2str(i+1),'}']; 
end
legend(legendInfo)

% ----------- Plot the original Trackmate data (Surface Plot) ------------
figure
surf(tint,xint,u','FaceAlpha',0.5)
xlabel('t (s)'); ylabel('x (\mu m)'); zlabel('u(x,t) (\mu m)');
title('TrackMate Surface Plot');
ax = gca;
ax.FontSize = 30; 
exportgraphics(ax,'TrackMateSurfacePlot.pdf','BackgroundColor','none')
view([-100 20])
ratio = (max(max(u))-min(min(u)))/(max(x_start)-min(x_start)); %setup for plotting since the scales aren't the same

%------------- Plot the original Trackmate data of u vs. x in time Animation -------------

% figure
% axis([min(x)-5*R max(x)+5*R min(min(X))-5*R max(max(X))+5*R])
% xlabel('x (\mu m)'); ylabel('u(x,t): Trackmate Displacement (\mu m)');
% title('Displacement Data from Trackmate - Animation');
% 
% %Plot the starting positions
% for i=1: length(x)
% PNA = [x(i),X(i,1)];
% c = [PNA(1),PNA(2)];
% pos1 = [c-R (1/ratio)*R 2*R];
% rect1 = rectangle('Position',pos1,'Curvature',[1,1], 'FaceColor','b');
% % pause(0.1)
% grid on
% end
% 
% % Plot the end points
% for i=1: length(x)
% PNA = [x(i),X(i,end)];
% c = [PNA(1),PNA(2)];
% pos1 = [c-R (1/ratio)*R 2*R];
% rect1 = rectangle('Position',pos1,'Curvature',[1,1], 'FaceColor','b');
% % pause(0.1)
% grid on
% end
% 
% hold on
% % Plot the trajectories over space and time
% for i=1: length(x)
% for j = 1: length(t)
% PNA = [x(i),X(i,j)];
% c = [PNA(1),PNA(2)];
% pos1 = [c-R (1/ratio)*R 2*R];
% rect1 = rectangle('Position',pos1,'Curvature',[1,1], 'FaceColor','b');
% pause(0.2)
% grid on
% delete(rect1)
% end
% end

% --------------------- Plot the original Trackmate data (Waterfall) ---------------------
u_trans = u';

figure
waterfall(tint,x_start(1:size(x_start,2)),u_trans(1:size(x_start,2),:))
xlabel('Time (seconds)'); ylabel('x (\mu m)'); zlabel('u(x,t): Trackmate Displacement (\mu m)');
title('Displacement Data from Trackmate');

% =============================================================================================================
%% Simulate the PDE of the System
mm = 0;
% xspace = linspace(-100,max(x_start),100);
xspace = linspace(min(x_start),max(x_start),100);
time = linspace(min(t),max(t),100);

sol = pdepe(mm,@pdex4pde,@pdex4ic,@pdex4bc,xspace,time);
ul = sol(:,:,1);
save('uPDE.mat','sol') 
save('uPDE.m','sol') 
save('uPDE.csv','sol') 
% --------------------- Plot the PDE Solution (Surface) ---------------------
% figure
hold on
xspace_selected = linspace(min(x_start),max(x_start),100);
surf(time,xspace,ul.','FaceAlpha',0.5)
ylabel('x (\mu m)'); xlabel('t (s)'); zlabel('u(x,t) (\mu m)');
title('Displacement from DMD')
ax = gca;
ax.FontSize = 30; 
exportgraphics(ax,'Project1_MasterFileSurfacePlot.pdf','BackgroundColor','none')
view([-100 20])
%--------------------------- Alternative Data Sets ------------------------------
% Original data set
% t_alt= [0;
% 1800.495546;
% 3600.991092;
% 5401.486637;
% 7201.982183;
% 9002.477729;
% 10802.97327;
% 12603.46882;
% 14403.96437;
% 16204.45991;
% 18004.95546;
% 19805.451;
% 21605.94655;
% 23406.4421;
% 25206.93764];
% 
% load x_cell.m;
% x_alt= sortrows(x_cell.');


% 1st alternative data set: x_cell_C2_R.m
% t_alt = [0;
% 1799.4076;
% 3598.815199;
% 5398.222799;
% 7197.630399;
% 8997.037999;
% 10796.4456;
% 12595.8532;
% 14395.2608;
% 16194.6684;
% 17994.076;
% 19793.4836;
% 21592.8912;
% 23392.2988;
% 25191.7064;
% 26991.114;
% 28790.5216;
% 30589.92919];
% % t_alt = 0:1799.4076:1799.4076*37;
% load x_cell_C2_R.m
% x_alt = sortrows(x_cell_C2_R.');

% 2nd alternative data set: x_cell_C2_R.m
t_alt = [0;
1799.4076;
3598.815199;
5398.222799;
8997.037999;
10796.4456;
12595.8532;
14395.2608;
16194.6684;
17994.076;
19793.4836;
21592.8912;
23392.2988;
25191.7064;
26991.114;
28790.5216;
30589.92919;
32389.33679;
34188.74439;
35988.15199;
37787.55959;
39586.96719;
41386.37479;
43185.78239;
44985.18999;
46784.59759;
48584.00519;
50383.41279;
52182.82039;
53982.22799;
55781.63559;
57581.04319;
59380.45079;
61179.85839;
62979.26599;
64778.67359;
66926.5379];
% t_alt = 0:1799.4076:1799.4076*37;
load x_cell_C2_L.m
x_alt = sortrows(x_cell_C2_L.');

for i = 1 : size(x_alt,2)
    for j = 1 : size(x_alt,1)
    u_alt(j,i) = x_alt(j,i) - x_alt(j,1);
    end
end

NGridPoints = 100;
x_start_alt = x_alt(:,1).'; %starting values of the cells
xspace_alt = linspace(min(x_start_alt),max(x_start_alt),NGridPoints);
time_alt = linspace(min(t_alt),max(t_alt),NGridPoints);

pos_lastcell_alt = size(u_alt,1);

% Extract the displacement of the RHS cell position
right_cell_polyfit_alt = polyfit(t_alt(1:end),u_alt(pos_lastcell_alt,:),PolyOrderOfApproxLastCell);
right_cell_polyval_alt = polyval(right_cell_polyfit_alt,time_alt(1:end));

% For C2_R (Data set 1)
% NumeratorR = u_alt(end,7)-u_alt(end,1);
% DenominatorR = (t_alt(7)-t_alt(1));
% InitialVelocity = NumeratorR/DenominatorR;
% DistanceToCoverCenterOfScratch = 216; 
% time_alt_Experimental = DistanceToCoverCenterOfScratch/(InitialVelocity);
% time_alt_Experimmental_vector = linspace(0,time_alt_Experimental,NGridPoints);

% For C2_L (Data set 2)
NumeratorR = u_alt(end,3)-u_alt(end,1);
DenominatorR = (t_alt(3)-t_alt(1));
InitialVelocity = NumeratorR/DenominatorR;
DistanceToCoverCenterOfScratch = 215; 
time_alt_Experimental = DistanceToCoverCenterOfScratch/(InitialVelocity);
time_alt_Experimmental_vector = linspace(0,time_alt_Experimental,NGridPoints);

mm = 0;
sol_alt = pdepe(mm,@pdex4pde,@pdex4ic,@pdex4bc_alt,xspace_alt,time_alt_Experimmental_vector);
ul_alt = sol_alt(:,:,1);

figure
waterfall(t_alt,x_start_alt(1:size(x_start_alt,2)),u_alt(1:size(x_start_alt,2),:))
xlabel('t (s)'); ylabel('x (\mu m)'); zlabel('u(x,t) (\mu m)');
title('Displacement from DMD'); % Displacement Data from Trackmate from Alternative Data Set #1
hold on
surf(time_alt_Experimmental_vector,xspace_alt,ul_alt.','FaceAlpha',0.5)
ax = gca;
ax.FontSize = 30; 
exportgraphics(ax,'Project1_MasterFileSurfacePlotDataSet2.pdf','BackgroundColor','none')
view([-100 20])

% figure()
% u_alt
% u_alt(ceil((x_start_alt(3)/(x_start_alt(end)-x_start_alt(1)))*length(xspace_alt)))

% error = 
% 
% ul_alt(50,:)
% u_alt(3,:)

% ul_alt(50,:)

abcd = linspace(x_start_alt(1),x_start_alt(end),NGridPoints)';
[c index] = min(abs(abcd-x_start_alt(1:size(x_start_alt,2))));
figure()
for i = 1 : size(u_alt,1)
    error(i,:) = u_alt(i,:)-interp1(time_alt,ul_alt(:,index(i))',linspace(1,time_alt(end),length(t_alt)));
    data(i) = mean(error(i,:));
    errhigh(i) = std(error(i,:));
    errlow(i) = std(error(i,:));
end
bar(x_start_alt,data)                
hold on
er = errorbar(x_start_alt,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
xlabel('Starting x-Value of the Cell'); ylabel('Average Error (\mu m)');title('Alternative Data Set #2');
ax = gca;
ax.FontSize = 34; 
exportgraphics(ax,'Project1_MasterFileStatistical3.pdf','BackgroundColor','none')

% sol in the workspace
u_94 = sol(:,94);

t = linspace(0,t(end),length(u_94));
dt = t(2);

%% Model Discovery
% State variables

figure(),plot(t,u_94,'Linewidth',[2])
xlabel('Time (seconds)'),ylabel('Displacement (\mu m): u_{94}(t) Obtained from TrackMate')
title('Displacement vs. Time')
% axis([0 25206.93764 -5 25])

% Compute the derivatives from this data
% This is the b in Ax=b
n = length(t);
for j = 2:n-1
    u_6dot(j-1) = (u_94(j+1)-u_94(j-1))/(2*dt);
end

% Note: We are only calculating the derivatives on the internal
%       points and will discard the external points. 
u_6s = u_94(2:n-1);
ts = t(2:n-1);

% A = [ones(length(u_6s),1),ts',(ts').*(ts'),(ts').*(ts').*(ts'),u_6s,u_6s.*u_6s,u_6s.*u_6s.*u_6s,u_6s.*u_6s.*u_6s.*u_6s];
A = [ones(length(u_6s),1),ts',(ts').*(ts'),(ts').*(ts').*(ts')];


%% Over-determined
% This system of equations can be easily solved by the \
xi1 = A\u_6dot.';

figure()
bar(xi1),xlabel('Terms in the A-matrix'),ylabel('$\dot{u}_{94}(t)$','Interpreter', 'latex')
% xticklabels({'1','t','t^{2}','t^{3}','u^{94}','(u_{94})^{2}','(u_{94})^{3}','(u_{94})^{4}'}) %From the A-matrix
xticklabels({'1','t','t^{2}','t^{3}'}) %From the A-matrix

% This system of equations can be easily solved by the pseudo-inverse
% xi1 = pinv(A)*u_6dot.';

% Uncomment the following in order to use the Lasso algorithm to promote sparsity
% xi1 = lasso(A,u_6dot.','Lambda',0.01);

%% ODE Solver for the SingleCellDyn
pos0 = 309.291;

dt =  0.01;
t = 0:dt:floor(t(end));

SV0 = [u(1)];

% Numerically integrate the displacement
[ts,SV] = ode45(@(ts,SV) SingleCellDyn(ts,SV,xi1),ts,SV0);
figure()
plot(ts,SV(:,1),'Linewidth',[2]') 
xlabel('Time (seconds)'),ylabel('Displacement (\mu m): u_{94}(t) Num. Int. from Discovered Dynamic Model')
title('Displacement vs. Time')

% Plot the position
figure()
plot(ts,SV(:,1) + pos0) 
xlabel('Time (seconds)'),ylabel('Position (\mu m): x_{94}')
title('Position vs. Time')

%% DMD of position
pos = SV(:,1) + pos0;
t = linspace(0,t(end),length(pos));
dt = t(2);
n = length(t);
for j = 2:n-1
    posdot(j-1) = (pos(j+1)-pos(j-1))/(2*dt);
end

poss = pos(2:n-1);
ts = t(2:n-1);

% B = [ones(length(poss),1),ts',(ts').*(ts'),(ts').*(ts').*(ts'),poss,poss.*poss,poss.*poss.*poss,poss.*poss.*poss.*poss];
B = [ones(length(poss),1),ts',(ts').*(ts'),(ts').*(ts').*(ts')];

xi2 = B\posdot.';

figure()
bar(xi2),xlabel('Terms in the A-matrix'),ylabel('$\dot{x}_{94}$','Interpreter', 'latex')
% xticklabels({'1','t','t^{2}','t^{3}','x_{94}','(x_{94})^{2}','(x_{94})^{3}','(x_{94})^{4}'}) %From the B-matrix
xticklabels({'1','t','t^{2}','t^{3}'}) %From the B-matrix

% Numerically integrate the displacement
[ts2,SV2] = ode45(@(ts,SV) SingleCellDyn2(ts,SV,xi2),ts,pos0);
figure()
plot(ts2,SV2(:,1)) 
xlabel('Time (seconds)'),ylabel('Position (\mu m): x_{94}')
title('Position vs. Time')

%% =============================== Functions ===============================
% ***
function rhs = SingleCellDyn(ts,SV,xi1)
% a = [1,ts,ts.*ts,ts.*ts.*ts,SV(1),SV(1).*SV(1),SV(1).*SV(1).*SV(1),SV(1).*SV(1).*SV(1).*SV(1)];
a = [1,ts,ts.*ts,ts.*ts.*ts];

rhs = (a*xi1)';
end

function rhs = SingleCellDyn2(ts,SV,xi2)
% b = [1,ts,ts.*ts,ts.*ts.*ts,SV(1),SV(1).*SV(1),SV(1).*SV(1).*SV(1),SV(1).*SV(1).*SV(1).*SV(1)];
b = [1,ts,ts.*ts,ts.*ts.*ts];

rhs = (b*xi2)';
end

function [c,f,s] = pdex4pde(x,t,u,DuDx)
global xi right_cell_polyfit right_cell_polyval left_cell_polyfit left_cell_polyval
c = 1;
f = xi(end)*DuDx;

% ***
% s=0;
% s = [u DuDx.*DuDx]*xi(1:end-1); %From the A-matrix 
% s = [DuDx,(DuDx)^2,(DuDx)^3,(DuDx)^4]*xi(7:10); %From the A-matrix 
% s = [(DuDx)^2]*xi(8); %From the A-matrix 

% s = [u u.^2 u.^3 u.^4 sin(u) sin(2*u) sin(u).^2 DuDx DuDx.*u DuDx.*DuDx]*xi(1:end-1); %From the A-matrix
% s = [u u.^2 u.^3 u.^4 sin(u) sin(2*u) sin(u).^2 DuDx DuDx.*u DuDx.*DuDx]*xi(1:end-1); %From the A-matrix

s = [DuDx]*xi(1:end-1);
% s = [DuDx (DuDx)^2]*xi(1:end-1);

% s = [1,u,u^2,u^3,u^4,u^5,DuDx,(DuDx)^2,(DuDx)^3,(DuDx)^4,DuDx*u,DuDx*u^2,DuDx*u^3,DuDx*u^4]*xi(1:end-1);
end

function u_0 = pdex4ic(x)
u_0=0;
end

function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,time)
global xi right_cell_polyfit right_cell_polyval left_cell_polyfit left_cell_polyval right_cell_polyfit_alt right_cell_polyval_alt

% has to take the format of the equation: p  + qf = 0

    for k = 1:size(right_cell_polyfit,2)
        timevec(k,1) = time^(k-1);
    end
    
% ***
%Left side: Dirichlet
% pl = ul -0;
% ql = 0;

% left side: Dirichlet
uL = left_cell_polyfit*flip(timevec); 
ql = 0;
pl = ul-uL;

% Left side: Neumann
% ql = 1;
% pl = 0;

%right side: Dirichlet (this one for the exact fit)
uR = right_cell_polyfit*flip(timevec);
qr=0;
pr = ur - uR;

% Right side: Dirichlet - (linearly interpolated to the mid of the scratch)
% uR = 39*time/(25206.93764);
% qr=0;
% pr = ur - uR;

end

function [pl,ql,pr,qr] = pdex4bc_alt(xl,ul,xr,ur,time)
global xi right_cell_polyfit right_cell_polyval left_cell_polyfit left_cell_polyval right_cell_polyfit_alt right_cell_polyval_alt
global u_alt t_alt
% has to take the format of the equation: p  + qf = 0

    for k = 1:size(right_cell_polyfit_alt,2)
        timevec(k,1) = time^(k-1);
    end
    
% ----------------- Dirichlet BC's Version #1 -----------------
% First Data Set
% NumeratorL = u_alt(1,7)-u_alt(1,1);
% DenominatorL = (t_alt(7)-t_alt(1));
% 
% uL  = (NumeratorL/DenominatorL)*time;
% ql=0;
% pl = ul - uL;   
% 
% NumeratorR = u_alt(end,7)-u_alt(end,1);
% DenominatorR = (t_alt(7)-t_alt(1));
% uR  = (NumeratorR/DenominatorR)*time;
% qr=0;
% pr = ur - uR;

% Second Data Set
NumeratorL = u_alt(1,3)-u_alt(1,1);
DenominatorL = (t_alt(3)-t_alt(1));

uL  = (NumeratorL/DenominatorL)*time;
ql=0;
pl = ul - uL;   

NumeratorR = u_alt(end,3)-u_alt(end,1);
DenominatorR = (t_alt(3)-t_alt(1));
uR  = (NumeratorR/DenominatorR)*time;
qr=0;
pr = ur - uR;
 
% ----------------- Dirichlet BC's Version #2 - Partial TrackMate Trajectories -----------------
% First Data Set:
% pl = ul - 0;
% ql = 0;
% 
% NumeratorR = u_alt(end,7)-u_alt(end,1);
% DenominatorR = (t_alt(7)-t_alt(1));
% uR  = (NumeratorR/DenominatorR)*time;
% qr=0;
% pr = ur - uR;

% Second Data Set:
% pl = ul - 0;
% ql = 0;
% 
% NumeratorR = u_alt(end,3)-u_alt(end,1);
% DenominatorR = (t_alt(3)-t_alt(1));
% uR  = (NumeratorR/DenominatorR)*time;
% qr=0;
% pr = ur - uR;

% ----------------- Mixed BCs (Dirichlet/Neumann) - Partial TrackMate Trajectories -----------------
% First Data Set:
% ql = 1;
% pl = 0;
% 
% NumeratorR = u_alt(end,7)-u_alt(end,1);
% DenominatorR = (t_alt(7)-t_alt(1));
% uR  = (NumeratorR/DenominatorR)*time;
% qr=0;
% pr = ur - uR;

% Second Data Set:
% ql = 1;
% pl = 0;
% 
% NumeratorR = u_alt(end,3)-u_alt(end,1);
% DenominatorR = (t_alt(3)-t_alt(1));
% uR  = (NumeratorR/DenominatorR)*time;
% qr=0;
% pr = ur - uR;
end
