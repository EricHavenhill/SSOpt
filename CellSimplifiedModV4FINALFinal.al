%----------------------------------------------------------------------
%    File: CellSimplifiedModV4FINALFinal.al
% Problem: Equations of Motion for the migrating cell
%          run CellSimplifiedModV4FINALFinal.al
%----------------------------------------------------------------------
% Coordinate system, bodies, lengths, etc.
Newtonian N
bodies A,B,C,D,E,F,G
constants L1 ,L2 ,L3 ,L4, L5, L6, L7, L8, L9
%----------------------------------------------------------------------
% Generalized coordinates
variables q{9}',u{9}'
q1' = u1
q2' = u2
q3' = u3
q4' = u4
q5' = u5
q6' = u6
q7' = u7
q8' = u8
q9' = u9
%----------------------------------------------------------------------
% Position and Orientation
simprot(N,A,3,q1) 
simprot(N,B,3,q1)

simprot(N,C,3,q3)
simprot(N,D,3,q3)

simprot(N,E,3,q5) 
simprot(N,F,3,q5) 

simprot(N,G,3,q7) 

P_No_Ao> = 0.5*L1*A1>
P_No_Bo> = q2*B1> + 0.5*L2*B1>

P_No_Co> = L9*N1> + 0.5*L3*C1>
P_No_Do> = L9*N1> + (q4+0.5*L4)*D1>

P_No_Eo> = 0.5*L9*N1> + 0.5*L5*E1>
P_No_Fo> = 0.5*L9*N1> + (q6+0.5*L6)*F1>

P_No_Go> = q8*N1> + q9*N2>
%----------------------------------------------------------------------
% Constraints

LOOP1> = q2*A1> + L2*A1> + 2*L7*G1> - (L4*D1>  + q4*D1> + L9*N1>)
LOOP2> = q2*A1> + L2*A1> + L7*G1> - L8*G2> - L6*F1> - q6*F1> - 0.5*L9*N1>
LOOP3> = q2*A1> + L2*A1> + L7*G1> - q9*N2> - q8*N1>

CONFIG[1]=dot(LOOP1>,N1>) 
CONFIG[2]=dot(LOOP1>,N2>) 
CONFIG[3]=dot(LOOP2>,N1>) 
CONFIG[4]=dot(LOOP2>,N2>) 
CONFIG[5]=dot(LOOP3>,N1>) 
CONFIG[6]=dot(LOOP3>,N2>) 

DEPENDENT=DT(CONFIG,N)
constrain(DEPENDENT[U1,U2,U3,U4,U5,U6]) 
%----------------------------------------------------------------------
% Velocities
W_A_N> = u1*A3> 
W_B_N> = u1*B3>
W_C_N> = u3*C3>
W_D_N> = u3*D3>
W_E_N> = u5*E3>
W_F_N> = u5*F3>
W_G_N> = u7*G3>

V_Ao_N> = dt(P_No_Ao>,N)
V_Bo_N> = dt(P_No_Bo>,N) 
V_Co_N> = dt(P_No_Co>,N) 
V_Do_N> = dt(P_No_Do>,N) 
V_Eo_N> = dt(P_No_Eo>,N) 
V_Fo_N> = dt(P_No_Fo>,N) 
V_Go_N> = explicit(dt(P_No_Go>,N),[u7, u8, u9])

%----------------------------------------------------------------------
% Accelerations

A_Ao_N> = dt(V_Ao_N>, N) 
A_Bo_N> = dt(V_Bo_N>, N) 
A_Co_N> = dt(V_Co_N>, N) 
A_Do_N> = dt(V_Do_N>, N) 
A_Eo_N> = dt(V_Eo_N>, N) 
A_Fo_N> = dt(V_Fo_N>, N) 
A_Go_N> = dt(V_Go_N>, N) 

ALF_A_N> = dt(W_A_N>, N) 
ALF_B_N> = dt(W_B_N>, N)
ALF_C_N> = dt(W_C_N>, N) 
ALF_D_N> = dt(W_D_N>, N) 
ALF_E_N> = dt(W_E_N>, N) 
ALF_F_N> = dt(W_F_N>, N) 
ALF_G_N> = dt(W_G_N>, N) 
%----------------------------------------------------------------------
% Mass & Inertia
mass A=mA
mass B=mB
mass C=mC
mass D=mD
mass E=mE
mass F=mF
mass G=mG

inertia A,IA1,IA2,IA3,0,0,0 
inertia B,IB1,IB2,IB3,0,0,0 
inertia C,IC1,IC2,IC3,0,0,0 
inertia D,ID1,ID2,ID3,0,0,0 
inertia E,IE1,IE2,IE3,0,0,0 
inertia F,IF1,IF2,IF3,0,0,0 
inertia G,IG1,IG2,IG3,0,0,0 
%----------------------------------------------------------------------
% Forces & Torques
constants grav 
gravity( -grav*N2> )

specified FB, FD, FF

Force_Bo> += FB*B1>
Force_Do> += FD*D1>
Force_Fo> += FF*F1>

ZERO = FR() + FRSTAR()

%----------------------------------------------------------------------
% EOM
% Note: The EOM are organized into the following form: 
%	A*q_dd + b + g = JTe*Fe + Gtran*Upsilon

A = -coef(zero,[u7',u8',u9'])
%A = [AA[2,2],AA[2,3];AA[3,2],AA[3,3]]

b = -exclude(zero,[u_independent';grav;FB;FD;FF]) 
%b = [bb[2];bb[3]]

g = -coef(zero,grav)*grav
%g = [gg[2];gg[3]]

GTran = coef(zero,[FB,FD,FF])
%GTranNew = [GTran[2,1],GTran[2,2];GTran[3,1],GTran[3,2]]

JT= coef(zero,[FB,FD,FF])
J=transpose(JT)
Jdot=DT(J)
detJ = det(J)

%jac=[dot(V_Go_N>,N1>);dot(V_Go_N>,N2>)]
%J=coef(jac,[u7,u8,u9])
%Jdot=dt(J)
%Jdotqdot=Jdot*[u7;u8;u9]
% ----------------------------------------------------------
% Control: Use the computed torque method to linearize the dynamics
UnitSystem kg, meter,sec
specified x8d,v8d,a8d
constants kv, kp, x7d,v7d,a7d,x9d,v9d,a9d

Input x7d = 0 rad, v7d = 0 rad/s, a7d = 0 rad/s^2
%Input x8d = 36.8624E-6 meter, v8d = 0 m/s, a8d = 0 m/s^2
Input x9d = 26E-6 meter, v9d = 0 m/s, a9d = 0 m/s^2

%x8d = 7.00E-6 + (7.12E-9)*T - ((1/2)*6.83E-13)*T^2 + ((1/3)*3.22E-17)*T^3 - ((1/4)*5.53E-22)*T^4
%v8d = (7.12E-9) - ((1/1)*1.37E-12)*T + ((1/1)*9.67E-17)*T^2 - ((1/1)*2.21E-21)*T^3
%a8d = -((1/1)*1.37E-12) + ((2/1)*9.67E-17)*T - ((3/1)*2.21E-21)*T^2

x8d = 7.00e-6 + 7.12e-9*t - 6.83e-13*t^2 + 3.22e-17*t^3 - 5.53e-22*t^4;
v8d = 7.12e-9 - 1.37e-12*t + 9.67e-17*t^2 - 2.21e-21*t^3;
a8d= -1.37e-12 + 1.93e-16*t - 6.63e-21*t^2;

UpStar = [a7d;a8d;a9d]+kv*([v7d;v8d;v9d]-[u7;u8;u9])+kp*([x7d;x8d;x9d]-[q7;dot(P_No_Go>,N1>);dot(P_No_Go>,N2>)])
Up = inverse(GTran)*(A*UpStar+b+g)
FB=Up[1]
FD=Up[2] 
FF=Up[3]
% ----------------------------------------------------------
% Input & Output
check = nicheck()

%Input tInitial=0, tFinal=1
Input tInitial=0, tFinal=25206.93764
%Input tInitial=0, tFinal=25200

% Input:
Input L1 = 20E-6 meter, L2 = 20E-6 meter, L3 = 20E-6 meter
Input L4 = 20E-6 meter, L5 = 20E-6 meter, L6 = 20E-6 meter
Input L7 = 7E-6 meter, L8 = 2.5E-6 meter, L9 = 60E-6 meter  

% Calc. the initial Q-values & place here
Input Q1 = 90 deg, Q2 = 6E-6 meter, Q3 = 150.58 deg
Input Q4 = 32.926E-6 meter, Q5 = 134.90 deg, Q6 = 12.759E-6 meter
Input Q7 = 0 rad, Q8 = 7E-6 meter, Q9 = 26E-6 meter

%Input U1 = 0 rad/s, U2 = 0 m/s, U3 = 0 rad/s
%Input U4 = 0 m/s, U5 = 0 rad/s, U6 = 0 m/s 
%Input U7 = 0 rad/s, U8 = 0 m/s, U9 = 0 m/s

Input MA = 1.23E-13 kg, MB = 1.23E-13 kg
Input MC = 1.23E-13 kg, MD = 1.23E-13 kg
Input ME = 1.23E-13 kg, MF = 1.23E-13 kg
Input MG = 2.29E-11 kg

Input grav = 9.81 m/sec^2

Input IA3 = (1/12)*(1.23E-13)*(20E-6)^2 kg*m^2
Input IB3 = (1/12)*(1.23E-13)*(20E-6)^2 kg*m^2
Input IC3 = (1/12)*(1.23E-13)*(20E-6)^2 kg*m^2
Input ID3 = (1/12)*(1.23E-13)*(20E-6)^2 kg*m^2
Input IE3 = (1/12)*(1.23E-13)*(20E-6)^2 kg*m^2
Input IF3 = (1/12)*(1.23E-13)*(20E-6)^2 kg*m^2
Input IG3 = (2/5)*(2.29E-11)*((7E-6)^2 +(2.5E-6)^2) kg*m^2

% Gains: kv=2*sqrt(kp)
Input kp = 1
Input kv = 2*sqrt(1)

Input integStp=1500

% Output to MATLAB
output t, q7,q8,q9,u7,u8,u9,FB,FD,FF,detJ
output A,b,g,GTran

save CellSimplifiedModV4FINALFinal.all
code Dynamics() CellSimplifiedModV4_Autolev_2_MATLAB.m