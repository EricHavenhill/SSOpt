function cellsimplifiedmodv4drugged_autolev_2_matlab
SolveOrdinaryDifferentialEquations
% File  cellsimplifiedmodv4drugged_autolev_2_matlab.m  created by Autolev 4.0 on Sat Oct 04 22:19:16 2025


%===========================================================================
function VAR = ReadUserInput
global   A7D A9D GRAV IA3 IB3 IC3 ID3 IE3 IF3 IG3 KP KV L1 L2 L3 L4 L5 L6 L7 L8 MA MB MC MD ME MF MG V7D V9D X7D X9D;
global   Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 U7 U8 U9 WCHECK1;
global   U1 U2 U3 U4 U5 U6 Q1p Q2p Q3p Q4p Q5p Q6p Q7p Q8p Q9p U7p U8p U9p WCHECK1p A8D FB FD FF V8D X8D;
global   DEGtoRAD RADtoDEG COEF RHS SolutionToAxEqualsB A B G GTRAN;
global   TINITIAL TFINAL INTEGSTP PRINTINT ABSERR RELERR;

%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
A7D                             =  0;                      % RAD/S^2             Constant
A9D                             =  0;                      % M/S^2               Constant
GRAV                            =  9.81;                   % M/SEC^2             Constant
IA3                             =  2.656666666666667E-26;  % KG*M^2              Constant
IB3                             =  2.656666666666667E-26;  % KG*M^2              Constant
IC3                             =  2.656666666666667E-26;  % KG*M^2              Constant
ID3                             =  2.656666666666667E-26;  % KG*M^2              Constant
IE3                             =  2.656666666666667E-26;  % KG*M^2              Constant
IF3                             =  2.656666666666667E-26;  % KG*M^2              Constant
IG3                             =  5.0609E-22;             % KG*M^2              Constant
KP                              =  1;                      % UNITS               Constant
KV                              =  2;                      % UNITS               Constant
L1                              =  20E-6;                  % METER               Constant
L2                              =  20E-6;                  % METER               Constant
L3                              =  20E-6;                  % METER               Constant
L4                              =  20E-6;                  % METER               Constant
L5                              =  20E-6;                  % METER               Constant
L6                              =  20E-6;                  % METER               Constant
L7                              =  7E-6;                   % METER               Constant
L8                              =  2.5E-6;                 % METER               Constant
MA                              =  7.97E-16;               % KG                  Constant
MB                              =  7.97E-16;               % KG                  Constant
MC                              =  7.97E-16;               % KG                  Constant
MD                              =  7.97E-16;               % KG                  Constant
ME                              =  7.97E-16;               % KG                  Constant
MF                              =  7.97E-16;               % KG                  Constant
MG                              =  2.29E-11;               % KG                  Constant
V7D                             =  0;                      % RAD/S               Constant
V9D                             =  0;                      % M/S                 Constant
X7D                             =  0;                      % RAD                 Constant
X9D                             =  26E-6;                  % METER               Constant

Q1                              =  90;                     % DEG                 Initial Value
Q2                              =  6E-6;                   % METER               Initial Value
Q3                              =  150.58;                 % DEG                 Initial Value
Q4                              =  32.926E-6;              % METER               Initial Value
Q5                              =  134.90;                 % DEG                 Initial Value
Q6                              =  12.759E-6;              % METER               Initial Value
Q7                              =  0;                      % RAD                 Initial Value
Q8                              =  7E-6;                   % METER               Initial Value
Q9                              =  26E-6;                  % METER               Initial Value
U7                              =  0.0;                    % UNITS               Initial Value
U8                              =  0.0;                    % UNITS               Initial Value
U9                              =  0.0;                    % UNITS               Initial Value
WCHECK1                         =  0.0;                    % UNITS               Initial Value

TINITIAL                        =  0;                      % SEC                 Initial Time
TFINAL                          =  23406.442;              % SEC                 Final Time
INTEGSTP                        =  1500;                   % SEC                 Integration Step
PRINTINT                        =  1;                      % Positive Integer    Print-Integer
ABSERR                          =  1.0E-08;                %                     Absolute Error
RELERR                          =  1.0E-07 ;               %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
Pi       = 3.141592653589793;
DEGtoRAD = Pi/180.0;
RADtoDEG = 180.0/Pi;
Q1 = Q1*DEGtoRAD;
Q3 = Q3*DEGtoRAD;
Q5 = Q5*DEGtoRAD;

% Reserve space and initialize matrices
COEF = zeros(3,3);
RHS = zeros(1,3);

% Evaluate constants
% Set the initial values of the states
VAR(1) = Q1;
VAR(2) = Q2;
VAR(3) = Q3;
VAR(4) = Q4;
VAR(5) = Q5;
VAR(6) = Q6;
VAR(7) = Q7;
VAR(8) = Q8;
VAR(9) = Q9;
VAR(10) = U7;
VAR(11) = U8;
VAR(12) = U9;
VAR(13) = WCHECK1;



%===========================================================================
function OpenOutputFilesAndWriteHeadings
FileIdentifier = fopen('cellsimplifiedmodv4drugged_autolev_2_matlab.1', 'wt');   if( FileIdentifier == -1 ) error('Error: unable to open file cellsimplifiedmodv4drugged_autolev_2_matlab.1'); end
fprintf( 1,             '%%       T             Q7             Q8             Q9             U7             U8             U9             FB             FD             FF           DET(J)\n' );
fprintf( 1,             '%%     (SEC)          (RAD)         (METER)        (METER)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
fprintf(FileIdentifier, '%% FILE: cellsimplifiedmodv4drugged_autolev_2_matlab.1\n%%\n' );
fprintf(FileIdentifier, '%%       T             Q7             Q8             Q9             U7             U8             U9             FB             FD             FF           DET(J)\n' );
fprintf(FileIdentifier, '%%     (SEC)          (RAD)         (METER)        (METER)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
FileIdentifier = fopen('cellsimplifiedmodv4drugged_autolev_2_matlab.2', 'wt');   if( FileIdentifier == -1 ) error('Error: unable to open file cellsimplifiedmodv4drugged_autolev_2_matlab.2'); end
fprintf(FileIdentifier, '%% FILE: cellsimplifiedmodv4drugged_autolev_2_matlab.2\n%%\n' );
fprintf(FileIdentifier, '%%\n' );
fprintf(FileIdentifier, '%%\n\n' );



%===========================================================================
% Main driver loop for numerical integration of differential equations
%===========================================================================
function SolveOrdinaryDifferentialEquations
global   A7D A9D GRAV IA3 IB3 IC3 ID3 IE3 IF3 IG3 KP KV L1 L2 L3 L4 L5 L6 L7 L8 MA MB MC MD ME MF MG V7D V9D X7D X9D;
global   Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 U7 U8 U9 WCHECK1;
global   U1 U2 U3 U4 U5 U6 Q1p Q2p Q3p Q4p Q5p Q6p Q7p Q8p Q9p U7p U8p U9p WCHECK1p A8D FB FD FF V8D X8D;
global   DEGtoRAD RADtoDEG COEF RHS SolutionToAxEqualsB A B G GTRAN;
global   TINITIAL TFINAL INTEGSTP PRINTINT ABSERR RELERR;

OpenOutputFilesAndWriteHeadings
VAR = ReadUserInput;
OdeMatlabOptions = odeset( 'RelTol',RELERR, 'AbsTol',ABSERR, 'MaxStep',INTEGSTP );
T = TINITIAL;
PrintCounter = 0;
mdlDerivatives(T,VAR,0);
while 1,
  if( TFINAL>=TINITIAL & T+0.01*INTEGSTP>=TFINAL ) PrintCounter = -1; end
  if( TFINAL<=TINITIAL & T+0.01*INTEGSTP<=TFINAL ) PrintCounter = -1; end
  if( PrintCounter <= 0.01 ),
     mdlOutputs(T,VAR,0);
     if( PrintCounter == -1 ) break; end
     PrintCounter = PRINTINT;
  end
  [TimeOdeArray,VarOdeArray] = ode45( @mdlDerivatives, [T T+INTEGSTP], VAR, OdeMatlabOptions, 0 );
  TimeAtEndOfArray = TimeOdeArray( length(TimeOdeArray) );
  if( abs(TimeAtEndOfArray - (T+INTEGSTP) ) >= abs(0.001*INTEGSTP) ) warning('numerical integration failed'); break;  end
  T = TimeAtEndOfArray;
  VAR = VarOdeArray( length(TimeOdeArray), : );
  PrintCounter = PrintCounter - 1;
end
mdlTerminate(T,VAR,0);



%===========================================================================
% mdlDerivatives: Calculates and returns the derivatives of the continuous states
%===========================================================================
function sys = mdlDerivatives(T,VAR,u)
global   A7D A9D GRAV IA3 IB3 IC3 ID3 IE3 IF3 IG3 KP KV L1 L2 L3 L4 L5 L6 L7 L8 MA MB MC MD ME MF MG V7D V9D X7D X9D;
global   Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 U7 U8 U9 WCHECK1;
global   U1 U2 U3 U4 U5 U6 Q1p Q2p Q3p Q4p Q5p Q6p Q7p Q8p Q9p U7p U8p U9p WCHECK1p A8D FB FD FF V8D X8D;
global   DEGtoRAD RADtoDEG COEF RHS SolutionToAxEqualsB A B G GTRAN;
global   TINITIAL TFINAL INTEGSTP PRINTINT ABSERR RELERR;

% Update variables after integration step
Q1 = VAR(1);
Q2 = VAR(2);
Q3 = VAR(3);
Q4 = VAR(4);
Q5 = VAR(5);
Q6 = VAR(6);
Q7 = VAR(7);
Q8 = VAR(8);
Q9 = VAR(9);
U7 = VAR(10);
U8 = VAR(11);
U9 = VAR(12);
WCHECK1 = VAR(13);
U1 = (cos(Q1)*U9-sin(Q1)*U8-L7*cos(Q1-Q7)*U7)/(L2+Q2);
Q1p = U1;
U2 = sin(Q1)*U9 + cos(Q1)*U8 - L7*sin(Q1-Q7)*U7;
Q2p = U2;
U3 = -(sin(Q3)*U8-cos(Q3)*U9-L7*cos(Q3-Q7)*U7)/(L4+Q4);
Q3p = U3;
U4 = sin(Q3)*U9 + cos(Q3)*U8 + L7*sin(Q3-Q7)*U7;
Q4p = U4;
U5 = (cos(Q5)*U9-sin(Q5)*U8-(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*U7)/(L6+Q6);
Q5p = U5;
U6 = sin(Q5)*U9 + cos(Q5)*U8 + (sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*U7;
Q6p = U6;
Q7p = U7;
Q8p = U8;
Q9p = U9;

% Quantities which were specified
A8D = 1.2E-12 + 8.28E-21*T^2 - 2.24E-16*T;
X8D = 7.0E-06 + 6.9E-22*T^4 + 6.0E-13*T^2 - 1.4E-09*T - 3.733333333333333E-17*T^3;
V8D = -1.4E-09 + 1.2E-12*T + 2.76E-21*T^3 - 1.12E-16*T^2;
FB = 0.25*((L7*cos(Q5)*sin(Q3-Q7)-cos(Q3)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))))*(2*MB*sin(Q1)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2)+2*MF*sin(Q5)*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2)+cos(Q1)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2)+cos(Q3)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4)+cos(Q5)*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6)-2*GRAV*(2*MG+2*MB*sin(Q1)^2+2*MD*sin(Q3)^2+2*MF*sin(Q5)^2+cos(Q1)^2*(L1*MA+MB*(L2+2*Q2))/(L2+Q2)+cos(Q3)^2*(L3*MC+MD*(L4+2*Q4))/(L4+Q4)+cos(Q5)^2*(L5*ME+MF*(L6+2*Q6))/(L6+Q6))-2*MD*sin(Q3)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2)-(4*MG+4*MB*sin(Q1)^2+4*MD*sin(Q3)^2+4*MF*sin(Q5)^2+cos(Q1)^2*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2+cos(Q3)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+cos(Q5)^2*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2)*(A9D+KP*(X9D-Q9)+KV*(V9D-U9))-(4*MB*sin(Q1)*cos(Q1)+4*MD*sin(Q3)*cos(Q3)+4*MF*sin(Q5)*cos(Q5)-sin(Q1)*cos(Q1)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-sin(Q3)*cos(Q3)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-sin(Q5)*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2)*(A8D+KP*(X8D-Q8)+KV*(V8D-U8))-(4*L7*MD*sin(Q3)*sin(Q3-Q7)+L7*cos(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+4*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-4*L7*MB*sin(Q1)*sin(Q1-Q7)-L7*cos(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A7D+KP*(X7D-Q7)+KV*(V7D-U7)))+(L7*sin(Q5)*sin(Q3-Q7)-sin(Q3)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))))*(2*MD*cos(Q3)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2)+sin(Q1)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2)+sin(Q3)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4)+sin(Q5)*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6)+(4*MB+4*MG+4*MD*cos(Q3)^2+4*MF*cos(Q5)^2+sin(Q3)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+sin(Q5)^2*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2-sin(Q1)^2*(4*MB-(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2))*(A8D+KP*(X8D-Q8)+KV*(V8D-U8))+(4*MB*sin(Q1)*cos(Q1)+4*MD*sin(Q3)*cos(Q3)+4*MF*sin(Q5)*cos(Q5)-sin(Q1)*cos(Q1)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-sin(Q3)*cos(Q3)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-sin(Q5)*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2)*(A9D+KP*(X9D-Q9)+KV*(V9D-U9))-2*GRAV*(L1*MA*sin(Q1)*cos(Q1)/(L2+Q2)+L3*MC*sin(Q3)*cos(Q3)/(L4+Q4)+L5*ME*sin(Q5)*cos(Q5)/(L6+Q6)-MB*sin(Q1)*cos(Q1)*(2-(L2+2*Q2)/(L2+Q2))-MD*sin(Q3)*cos(Q3)*(2-(L4+2*Q4)/(L4+Q4))-MF*sin(Q5)*cos(Q5)*(2-(L6+2*Q6)/(L6+Q6)))-2*MB*cos(Q1)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2)-2*MF*cos(Q5)*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2)-(4*L7*MB*cos(Q1)*sin(Q1-Q7)+L7*sin(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-4*L7*MD*cos(Q3)*sin(Q3-Q7)-L7*sin(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-4*MF*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-sin(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A7D+KP*(X7D-Q7)+KV*(V7D-U7)))-sin(Q3-Q5)*(2*MF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2)+L7*cos(Q3-Q7)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4)+(4*L7*MB*cos(Q1)*sin(Q1-Q7)+L7*sin(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-4*L7*MD*cos(Q3)*sin(Q3-Q7)-L7*sin(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-4*MF*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-sin(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A8D+KP*(X8D-Q8)+KV*(V8D-U8))-2*GRAV*(L3*L7*MC*cos(Q3)*cos(Q3-Q7)/(L4+Q4)+L7*MD*(2*sin(Q3)*sin(Q3-Q7)+cos(Q3)*(L4+2*Q4)*cos(Q3-Q7)/(L4+Q4))+2*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-L1*L7*MA*cos(Q1)*cos(Q1-Q7)/(L2+Q2)-L7*MB*(2*sin(Q1)*sin(Q1-Q7)+cos(Q1)*(L2+2*Q2)*cos(Q1-Q7)/(L2+Q2))-cos(Q5)*(L5*ME+MF*(L6+2*Q6))*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6))-2*L7*MB*sin(Q1-Q7)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2)-2*L7*MD*sin(Q3-Q7)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2)-L7*cos(Q1-Q7)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2)-(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6)-(4*IG3+4*MB*L7^2*sin(Q1-Q7)^2+4*MD*L7^2*sin(Q3-Q7)^2+L7^2*cos(Q1-Q7)^2*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2+L7^2*cos(Q3-Q7)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+4*MF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2+(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2/(L6+Q6)^2)*(A7D+KP*(X7D-Q7)+KV*(V7D-U7))-(4*L7*MD*sin(Q3)*sin(Q3-Q7)+L7*cos(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+4*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-4*L7*MB*sin(Q1)*sin(Q1-Q7)-L7*cos(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A9D+KP*(X9D-Q9)+KV*(V9D-U9))))/(L7*sin(Q5)*(cos(Q1)*sin(Q3-Q7)+cos(Q3)*sin(Q1-Q7))-sin(Q3)*(L7*cos(Q5)*sin(Q1-Q7)+cos(Q1)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))))-sin(Q1)*(L7*cos(Q5)*sin(Q3-Q7)-cos(Q3)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))));
FD = 0.25*((L7*cos(Q5)*sin(Q1-Q7)+cos(Q1)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))))*(2*MB*sin(Q1)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2)+2*MF*sin(Q5)*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2)+cos(Q1)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2)+cos(Q3)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4)+cos(Q5)*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6)-2*GRAV*(2*MG+2*MB*sin(Q1)^2+2*MD*sin(Q3)^2+2*MF*sin(Q5)^2+cos(Q1)^2*(L1*MA+MB*(L2+2*Q2))/(L2+Q2)+cos(Q3)^2*(L3*MC+MD*(L4+2*Q4))/(L4+Q4)+cos(Q5)^2*(L5*ME+MF*(L6+2*Q6))/(L6+Q6))-2*MD*sin(Q3)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2)-(4*MG+4*MB*sin(Q1)^2+4*MD*sin(Q3)^2+4*MF*sin(Q5)^2+cos(Q1)^2*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2+cos(Q3)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+cos(Q5)^2*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2)*(A9D+KP*(X9D-Q9)+KV*(V9D-U9))-(4*MB*sin(Q1)*cos(Q1)+4*MD*sin(Q3)*cos(Q3)+4*MF*sin(Q5)*cos(Q5)-sin(Q1)*cos(Q1)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-sin(Q3)*cos(Q3)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-sin(Q5)*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2)*(A8D+KP*(X8D-Q8)+KV*(V8D-U8))-(4*L7*MD*sin(Q3)*sin(Q3-Q7)+L7*cos(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+4*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-4*L7*MB*sin(Q1)*sin(Q1-Q7)-L7*cos(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A7D+KP*(X7D-Q7)+KV*(V7D-U7)))+(L7*sin(Q5)*sin(Q1-Q7)+sin(Q1)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))))*(2*MD*cos(Q3)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2)+sin(Q1)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2)+sin(Q3)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4)+sin(Q5)*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6)+(4*MB+4*MG+4*MD*cos(Q3)^2+4*MF*cos(Q5)^2+sin(Q3)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+sin(Q5)^2*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2-sin(Q1)^2*(4*MB-(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2))*(A8D+KP*(X8D-Q8)+KV*(V8D-U8))+(4*MB*sin(Q1)*cos(Q1)+4*MD*sin(Q3)*cos(Q3)+4*MF*sin(Q5)*cos(Q5)-sin(Q1)*cos(Q1)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-sin(Q3)*cos(Q3)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-sin(Q5)*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2)*(A9D+KP*(X9D-Q9)+KV*(V9D-U9))-2*GRAV*(L1*MA*sin(Q1)*cos(Q1)/(L2+Q2)+L3*MC*sin(Q3)*cos(Q3)/(L4+Q4)+L5*ME*sin(Q5)*cos(Q5)/(L6+Q6)-MB*sin(Q1)*cos(Q1)*(2-(L2+2*Q2)/(L2+Q2))-MD*sin(Q3)*cos(Q3)*(2-(L4+2*Q4)/(L4+Q4))-MF*sin(Q5)*cos(Q5)*(2-(L6+2*Q6)/(L6+Q6)))-2*MB*cos(Q1)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2)-2*MF*cos(Q5)*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2)-(4*L7*MB*cos(Q1)*sin(Q1-Q7)+L7*sin(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-4*L7*MD*cos(Q3)*sin(Q3-Q7)-L7*sin(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-4*MF*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-sin(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A7D+KP*(X7D-Q7)+KV*(V7D-U7)))+sin(Q1-Q5)*(2*MF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2)+L7*cos(Q3-Q7)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4)+(4*L7*MB*cos(Q1)*sin(Q1-Q7)+L7*sin(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-4*L7*MD*cos(Q3)*sin(Q3-Q7)-L7*sin(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-4*MF*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-sin(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A8D+KP*(X8D-Q8)+KV*(V8D-U8))-2*GRAV*(L3*L7*MC*cos(Q3)*cos(Q3-Q7)/(L4+Q4)+L7*MD*(2*sin(Q3)*sin(Q3-Q7)+cos(Q3)*(L4+2*Q4)*cos(Q3-Q7)/(L4+Q4))+2*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-L1*L7*MA*cos(Q1)*cos(Q1-Q7)/(L2+Q2)-L7*MB*(2*sin(Q1)*sin(Q1-Q7)+cos(Q1)*(L2+2*Q2)*cos(Q1-Q7)/(L2+Q2))-cos(Q5)*(L5*ME+MF*(L6+2*Q6))*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6))-2*L7*MB*sin(Q1-Q7)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2)-2*L7*MD*sin(Q3-Q7)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2)-L7*cos(Q1-Q7)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2)-(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6)-(4*IG3+4*MB*L7^2*sin(Q1-Q7)^2+4*MD*L7^2*sin(Q3-Q7)^2+L7^2*cos(Q1-Q7)^2*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2+L7^2*cos(Q3-Q7)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+4*MF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2+(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2/(L6+Q6)^2)*(A7D+KP*(X7D-Q7)+KV*(V7D-U7))-(4*L7*MD*sin(Q3)*sin(Q3-Q7)+L7*cos(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+4*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-4*L7*MB*sin(Q1)*sin(Q1-Q7)-L7*cos(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A9D+KP*(X9D-Q9)+KV*(V9D-U9))))/(L7*sin(Q5)*(cos(Q1)*sin(Q3-Q7)+cos(Q3)*sin(Q1-Q7))-sin(Q3)*(L7*cos(Q5)*sin(Q1-Q7)+cos(Q1)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))))-sin(Q1)*(L7*cos(Q5)*sin(Q3-Q7)-cos(Q3)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))));
FF = -0.25*(L7*(cos(Q1)*sin(Q3-Q7)+cos(Q3)*sin(Q1-Q7))*(2*MB*sin(Q1)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2)+2*MF*sin(Q5)*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2)+cos(Q1)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2)+cos(Q3)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4)+cos(Q5)*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6)-2*GRAV*(2*MG+2*MB*sin(Q1)^2+2*MD*sin(Q3)^2+2*MF*sin(Q5)^2+cos(Q1)^2*(L1*MA+MB*(L2+2*Q2))/(L2+Q2)+cos(Q3)^2*(L3*MC+MD*(L4+2*Q4))/(L4+Q4)+cos(Q5)^2*(L5*ME+MF*(L6+2*Q6))/(L6+Q6))-2*MD*sin(Q3)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2)-(4*MG+4*MB*sin(Q1)^2+4*MD*sin(Q3)^2+4*MF*sin(Q5)^2+cos(Q1)^2*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2+cos(Q3)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+cos(Q5)^2*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2)*(A9D+KP*(X9D-Q9)+KV*(V9D-U9))-(4*MB*sin(Q1)*cos(Q1)+4*MD*sin(Q3)*cos(Q3)+4*MF*sin(Q5)*cos(Q5)-sin(Q1)*cos(Q1)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-sin(Q3)*cos(Q3)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-sin(Q5)*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2)*(A8D+KP*(X8D-Q8)+KV*(V8D-U8))-(4*L7*MD*sin(Q3)*sin(Q3-Q7)+L7*cos(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+4*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-4*L7*MB*sin(Q1)*sin(Q1-Q7)-L7*cos(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A7D+KP*(X7D-Q7)+KV*(V7D-U7)))+L7*(sin(Q1)*sin(Q3-Q7)+sin(Q3)*sin(Q1-Q7))*(2*MD*cos(Q3)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2)+sin(Q1)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2)+sin(Q3)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4)+sin(Q5)*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6)+(4*MB+4*MG+4*MD*cos(Q3)^2+4*MF*cos(Q5)^2+sin(Q3)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+sin(Q5)^2*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2-sin(Q1)^2*(4*MB-(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2))*(A8D+KP*(X8D-Q8)+KV*(V8D-U8))+(4*MB*sin(Q1)*cos(Q1)+4*MD*sin(Q3)*cos(Q3)+4*MF*sin(Q5)*cos(Q5)-sin(Q1)*cos(Q1)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-sin(Q3)*cos(Q3)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-sin(Q5)*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2)*(A9D+KP*(X9D-Q9)+KV*(V9D-U9))-2*GRAV*(L1*MA*sin(Q1)*cos(Q1)/(L2+Q2)+L3*MC*sin(Q3)*cos(Q3)/(L4+Q4)+L5*ME*sin(Q5)*cos(Q5)/(L6+Q6)-MB*sin(Q1)*cos(Q1)*(2-(L2+2*Q2)/(L2+Q2))-MD*sin(Q3)*cos(Q3)*(2-(L4+2*Q4)/(L4+Q4))-MF*sin(Q5)*cos(Q5)*(2-(L6+2*Q6)/(L6+Q6)))-2*MB*cos(Q1)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2)-2*MF*cos(Q5)*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2)-(4*L7*MB*cos(Q1)*sin(Q1-Q7)+L7*sin(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-4*L7*MD*cos(Q3)*sin(Q3-Q7)-L7*sin(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-4*MF*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-sin(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A7D+KP*(X7D-Q7)+KV*(V7D-U7)))+sin(Q1-Q3)*(2*MF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2)+L7*cos(Q3-Q7)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4)+(4*L7*MB*cos(Q1)*sin(Q1-Q7)+L7*sin(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2-4*L7*MD*cos(Q3)*sin(Q3-Q7)-L7*sin(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-4*MF*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-sin(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A8D+KP*(X8D-Q8)+KV*(V8D-U8))-2*GRAV*(L3*L7*MC*cos(Q3)*cos(Q3-Q7)/(L4+Q4)+L7*MD*(2*sin(Q3)*sin(Q3-Q7)+cos(Q3)*(L4+2*Q4)*cos(Q3-Q7)/(L4+Q4))+2*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-L1*L7*MA*cos(Q1)*cos(Q1-Q7)/(L2+Q2)-L7*MB*(2*sin(Q1)*sin(Q1-Q7)+cos(Q1)*(L2+2*Q2)*cos(Q1-Q7)/(L2+Q2))-cos(Q5)*(L5*ME+MF*(L6+2*Q6))*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6))-2*L7*MB*sin(Q1-Q7)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2)-2*L7*MD*sin(Q3-Q7)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2)-L7*cos(Q1-Q7)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2)-(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6)-(4*IG3+4*MB*L7^2*sin(Q1-Q7)^2+4*MD*L7^2*sin(Q3-Q7)^2+L7^2*cos(Q1-Q7)^2*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2+L7^2*cos(Q3-Q7)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+4*MF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2+(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2/(L6+Q6)^2)*(A7D+KP*(X7D-Q7)+KV*(V7D-U7))-(4*L7*MD*sin(Q3)*sin(Q3-Q7)+L7*cos(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2+4*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-4*L7*MB*sin(Q1)*sin(Q1-Q7)-L7*cos(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2-cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2)*(A9D+KP*(X9D-Q9)+KV*(V9D-U9))))/(L7*sin(Q5)*(cos(Q1)*sin(Q3-Q7)+cos(Q3)*sin(Q1-Q7))-sin(Q3)*(L7*cos(Q5)*sin(Q1-Q7)+cos(Q1)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))))-sin(Q1)*(L7*cos(Q5)*sin(Q3-Q7)-cos(Q3)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))));

WCHECK1p = 0.5*(GRAV*L3*L7*MC*cos(Q3)*cos(Q3-Q7)/(L4+Q4)+2*GRAV*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))+L7*(2*FB*sin(Q1-Q7)-2*GRAV*MB*sin(Q1)*sin(Q1-Q7)-GRAV*MB*cos(Q1)*(L2+2*Q2)*cos(Q1-Q7)/(L2+Q2))-GRAV*L1*L7*MA*cos(Q1)*cos(Q1-Q7)/(L2+Q2)-2*FF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-L7*(2*FD*sin(Q3-Q7)-2*GRAV*MD*sin(Q3)*sin(Q3-Q7)-GRAV*MD*cos(Q3)*(L4+2*Q4)*cos(Q3-Q7)/(L4+Q4))-GRAV*cos(Q5)*(L5*ME+MF*(L6+2*Q6))*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6))*U7 - 0.5*(2*FD*sin(Q3)+2*FF*sin(Q5)-2*GRAV*MG-2*GRAV*MD*sin(Q3)^2-2*GRAV*MF*sin(Q5)^2-GRAV*(L1*MA+MB*(L2+2*Q2))/(L2+Q2)-GRAV*cos(Q3)^2*(L3*MC+MD*(L4+2*Q4))/(L4+Q4)-GRAV*cos(Q5)^2*(L5*ME+MF*(L6+2*Q6))/(L6+Q6)-sin(Q1)*(2*GRAV*MB*sin(Q1)-2*FB-GRAV*sin(Q1)*(L1*MA+MB*(L2+2*Q2))/(L2+Q2)))*U9 - 0.5*(GRAV*L1*MA*sin(Q1)*cos(Q1)/(L2+Q2)+GRAV*L3*MC*sin(Q3)*cos(Q3)/(L4+Q4)+GRAV*L5*ME*sin(Q5)*cos(Q5)/(L6+Q6)-cos(Q1)*(2*GRAV*MB*sin(Q1)-2*FB-GRAV*MB*sin(Q1)*(L2+2*Q2)/(L2+Q2))-cos(Q3)*(2*GRAV*MD*sin(Q3)-2*FD-GRAV*MD*sin(Q3)*(L4+2*Q4)/(L4+Q4))-cos(Q5)*(2*GRAV*MF*sin(Q5)-2*FF-GRAV*MF*sin(Q5)*(L6+2*Q6)/(L6+Q6)))*U8;

COEF(1,1) = -IG3 - MB*L7^2*sin(Q1-Q7)^2 - MD*L7^2*sin(Q3-Q7)^2 - 0.25*L7^2*cos(Q1-Q7)^2*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 - 0.25*L7^2*cos(Q3-Q7)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 - MF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2 - 0.25*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2/(L6+Q6)^2;
COEF(1,2) = L7*MB*cos(Q1)*sin(Q1-Q7) + 0.25*L7*sin(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 - L7*MD*cos(Q3)*sin(Q3-Q7) - 0.25*L7*sin(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 - MF*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) - 0.25*sin(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2;
COEF(1,3) = L7*MB*sin(Q1)*sin(Q1-Q7) + 0.25*L7*cos(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 + 0.25*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2 - L7*MD*sin(Q3)*sin(Q3-Q7) - 0.25*L7*cos(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 - MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)));
COEF(2,1) = L7*MB*cos(Q1)*sin(Q1-Q7) + 0.25*L7*sin(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 - L7*MD*cos(Q3)*sin(Q3-Q7) - 0.25*L7*sin(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 - MF*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) - 0.25*sin(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2;
COEF(2,2) = 0.25*sin(Q1)^2*(4*MB-(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2) - MB - MG - MD*cos(Q3)^2 - MF*cos(Q5)^2 - 0.25*sin(Q3)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 - 0.25*sin(Q5)^2*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2;
COEF(2,3) = 0.25*sin(Q1)*cos(Q1)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 + 0.25*sin(Q3)*cos(Q3)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 + 0.25*sin(Q5)*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2 - MB*sin(Q1)*cos(Q1) - MD*sin(Q3)*cos(Q3) - MF*sin(Q5)*cos(Q5);
COEF(3,1) = L7*MB*sin(Q1)*sin(Q1-Q7) + 0.25*L7*cos(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 + 0.25*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2 - L7*MD*sin(Q3)*sin(Q3-Q7) - 0.25*L7*cos(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 - MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)));
COEF(3,2) = 0.25*sin(Q1)*cos(Q1)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 + 0.25*sin(Q3)*cos(Q3)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 + 0.25*sin(Q5)*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2 - MB*sin(Q1)*cos(Q1) - MD*sin(Q3)*cos(Q3) - MF*sin(Q5)*cos(Q5);
COEF(3,3) = -MG - MB*sin(Q1)^2 - MD*sin(Q3)^2 - MF*sin(Q5)^2 - 0.25*cos(Q1)^2*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 - 0.25*cos(Q3)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 - 0.25*cos(Q5)^2*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2;
RHS(1) = GRAV*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) + 0.5*L7*(2*FB*sin(Q1-Q7)-2*GRAV*MB*sin(Q1)*sin(Q1-Q7)-GRAV*MB*cos(Q1)*(L2+2*Q2)*cos(Q1-Q7)/(L2+Q2)) + 0.5*L7*MB*sin(Q1-Q7)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2) + 0.5*L7*MD*sin(Q3-Q7)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2) + 0.25*L7*cos(Q1-Q7)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-2*GRAV*L1*MA*cos(Q1)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2) + 0.25*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-2*GRAV*L5*ME*cos(Q5)-2*GRAV*MF*cos(Q5)*(L6+2*Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6) - FF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) - 0.5*L7*(2*FD*sin(Q3-Q7)-2*GRAV*MD*sin(Q3)*sin(Q3-Q7)-GRAV*MD*cos(Q3)*(L4+2*Q4)*cos(Q3-Q7)/(L4+Q4)) - 0.5*MF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2) - 0.25*L7*cos(Q3-Q7)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-2*GRAV*L3*MC*cos(Q3)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4);
RHS(2) = 0.5*cos(Q1)*(2*GRAV*MB*sin(Q1)-2*FB-GRAV*MB*sin(Q1)*(L2+2*Q2)/(L2+Q2)) + 0.5*cos(Q3)*(2*GRAV*MD*sin(Q3)-2*FD-GRAV*MD*sin(Q3)*(L4+2*Q4)/(L4+Q4)) + 0.5*cos(Q5)*(2*GRAV*MF*sin(Q5)-2*FF-GRAV*MF*sin(Q5)*(L6+2*Q6)/(L6+Q6)) + 0.5*MD*cos(Q3)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2) + 0.25*sin(Q1)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-2*GRAV*L1*MA*cos(Q1)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2) + 0.25*sin(Q3)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-2*GRAV*L3*MC*cos(Q3)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4) + 0.25*sin(Q5)*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-2*GRAV*L5*ME*cos(Q5)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6) - 0.5*MB*cos(Q1)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2) - 0.5*MF*cos(Q5)*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2);
RHS(3) = GRAV*MG + GRAV*MB*sin(Q1)^2 + GRAV*MD*sin(Q3)^2 + GRAV*MF*sin(Q5)^2 + 0.5*MD*sin(Q3)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2) - FB*sin(Q1) - FD*sin(Q3) - FF*sin(Q5) - 0.5*MB*sin(Q1)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2) - 0.5*MF*sin(Q5)*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2) - 0.25*cos(Q1)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-2*GRAV*L1*MA*cos(Q1)-2*GRAV*MB*cos(Q1)*(L2+2*Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2) - 0.25*cos(Q3)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-2*GRAV*L3*MC*cos(Q3)-2*GRAV*MD*cos(Q3)*(L4+2*Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4) - 0.25*cos(Q5)*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-2*GRAV*L5*ME*cos(Q5)-2*GRAV*MF*cos(Q5)*(L6+2*Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6);
SolutionToAxEqualsB = COEF\RHS';

% Update variables after uncoupling equations
U7p = SolutionToAxEqualsB(1);
U8p = SolutionToAxEqualsB(2);
U9p = SolutionToAxEqualsB(3);

% Update derivative array prior to integration step
VARp(1) = Q1p;
VARp(2) = Q2p;
VARp(3) = Q3p;
VARp(4) = Q4p;
VARp(5) = Q5p;
VARp(6) = Q6p;
VARp(7) = Q7p;
VARp(8) = Q8p;
VARp(9) = Q9p;
VARp(10) = U7p;
VARp(11) = U8p;
VARp(12) = U9p;
VARp(13) = WCHECK1p;

sys = VARp';



%===========================================================================
% mdlOutputs: Calculates and return the outputs
%===========================================================================
function Output = mdlOutputs(T,VAR,u)
global   A7D A9D GRAV IA3 IB3 IC3 ID3 IE3 IF3 IG3 KP KV L1 L2 L3 L4 L5 L6 L7 L8 MA MB MC MD ME MF MG V7D V9D X7D X9D;
global   Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 U7 U8 U9 WCHECK1;
global   U1 U2 U3 U4 U5 U6 Q1p Q2p Q3p Q4p Q5p Q6p Q7p Q8p Q9p U7p U8p U9p WCHECK1p A8D FB FD FF V8D X8D;
global   DEGtoRAD RADtoDEG COEF RHS SolutionToAxEqualsB A B G GTRAN;
global   TINITIAL TFINAL INTEGSTP PRINTINT ABSERR RELERR;

% Evaluate output quantities
Output(1)=T;  Output(2)=Q7;  Output(3)=Q8;  Output(4)=Q9;  Output(5)=U7;  Output(6)=U8;  Output(7)=U9;  Output(8)=FB;  Output(9)=FD;  Output(10)=FF;  Output(11)=(L7*cos(Q5)*(sin(Q1)*sin(Q3-Q7)+sin(Q3)*sin(Q1-Q7))-L7*sin(Q5)*(cos(Q1)*sin(Q3-Q7)+cos(Q3)*sin(Q1-Q7))-sin(Q1-Q3)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))));
Output(12)=0.0;
  Output(13)=0.0;
  Output(14)=0.0;
  Output(15)=0.0;

FileIdentifier = fopen('all');
WriteOutput( 1,                 Output(1:11) );
WriteOutput( FileIdentifier(1), Output(1:11) );
WriteOutput( FileIdentifier(2), Output(12:15) );
A(1,1) = IG3 + MB*L7^2*sin(Q1-Q7)^2 + MD*L7^2*sin(Q3-Q7)^2 + 0.25*L7^2*cos(Q1-Q7)^2*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 + 0.25*L7^2*cos(Q3-Q7)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 + MF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2 + 0.25*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2/(L6+Q6)^2;
A(1,2) = L7*MD*cos(Q3)*sin(Q3-Q7) + 0.25*L7*sin(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 + MF*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) + 0.25*sin(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2 - L7*MB*cos(Q1)*sin(Q1-Q7) - 0.25*L7*sin(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2;
A(1,3) = L7*MD*sin(Q3)*sin(Q3-Q7) + 0.25*L7*cos(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 + MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) - L7*MB*sin(Q1)*sin(Q1-Q7) - 0.25*L7*cos(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 - 0.25*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2;
A(2,1) = L7*MD*cos(Q3)*sin(Q3-Q7) + 0.25*L7*sin(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 + MF*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) + 0.25*sin(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2 - L7*MB*cos(Q1)*sin(Q1-Q7) - 0.25*L7*sin(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2;
A(2,2) = MB + MG + MD*cos(Q3)^2 + MF*cos(Q5)^2 + 0.25*sin(Q3)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 + 0.25*sin(Q5)^2*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2 - 0.25*sin(Q1)^2*(4*MB-(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2);
A(2,3) = MB*sin(Q1)*cos(Q1) + MD*sin(Q3)*cos(Q3) + MF*sin(Q5)*cos(Q5) - 0.25*sin(Q1)*cos(Q1)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 - 0.25*sin(Q3)*cos(Q3)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 - 0.25*sin(Q5)*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2;
A(3,1) = L7*MD*sin(Q3)*sin(Q3-Q7) + 0.25*L7*cos(Q3)*cos(Q3-Q7)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 + MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) - L7*MB*sin(Q1)*sin(Q1-Q7) - 0.25*L7*cos(Q1)*cos(Q1-Q7)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 - 0.25*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6)^2;
A(3,2) = MB*sin(Q1)*cos(Q1) + MD*sin(Q3)*cos(Q3) + MF*sin(Q5)*cos(Q5) - 0.25*sin(Q1)*cos(Q1)*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 - 0.25*sin(Q3)*cos(Q3)*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 - 0.25*sin(Q5)*cos(Q5)*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2;
A(3,3) = MG + MB*sin(Q1)^2 + MD*sin(Q3)^2 + MF*sin(Q5)^2 + 0.25*cos(Q1)^2*(4*IA3+4*IB3+MA*L1^2+MB*(L2+2*Q2)^2)/(L2+Q2)^2 + 0.25*cos(Q3)^2*(4*IC3+4*ID3+MC*L3^2+MD*(L4+2*Q4)^2)/(L4+Q4)^2 + 0.25*cos(Q5)^2*(4*IE3+4*IF3+ME*L5^2+MF*(L6+2*Q6)^2)/(L6+Q6)^2;
B(1) = 0.5*L7*MB*sin(Q1-Q7)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2) + 0.5*L7*MD*sin(Q3-Q7)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2) + 0.25*L7*cos(Q1-Q7)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2) + 0.25*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6) - 0.5*MF*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2) - 0.25*L7*cos(Q3-Q7)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4);
B(2) = 0.5*MD*cos(Q3)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2) + 0.25*sin(Q1)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2) + 0.25*sin(Q3)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4) + 0.25*sin(Q5)*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6) - 0.5*MB*cos(Q1)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2) - 0.5*MF*cos(Q5)*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2);
B(3) = 0.5*MD*sin(Q3)*(2*(L4+Q4)*U3^2-(L4+2*Q4)*U3^2-2*L7*cos(Q3-Q7)*U7^2) - 0.5*MB*sin(Q1)*((L2+2*Q2)*U1^2-2*(L2+Q2)*U1^2-2*L7*cos(Q1-Q7)*U7^2) - 0.5*MF*sin(Q5)*((L6+2*Q6)*U5^2-2*(L6+Q6)*U5^2-2*L8*sin(Q5-Q7)*U7^2) - 0.25*cos(Q1)*((4*IA3+4*IB3+MA*L1^2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)-MB*(L2+2*Q2)*(4*U1*U2-(L2+2*Q2)*(2*U1*U2+L7*sin(Q1-Q7)*U7^2)/(L2+Q2)))/(L2+Q2) - 0.25*cos(Q3)*((4*IC3+4*ID3+MC*L3^2)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)-MD*(L4+2*Q4)*(4*U3*U4-(L4+2*Q4)*(2*U3*U4-L7*sin(Q3-Q7)*U7^2)/(L4+Q4)))/(L4+Q4) - 0.25*cos(Q5)*((4*IE3+4*IF3+ME*L5^2)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)-MF*(L6+2*Q6)*(4*U5*U6-(L6+2*Q6)*(2*U5*U6-L8*cos(Q5-Q7)*U7^2)/(L6+Q6)))/(L6+Q6);
G(1) = 0.5*GRAV*(L3*L7*MC*cos(Q3)*cos(Q3-Q7)/(L4+Q4)+L7*MD*(2*sin(Q3)*sin(Q3-Q7)+cos(Q3)*(L4+2*Q4)*cos(Q3-Q7)/(L4+Q4))+2*MF*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin(Q5-Q7)-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-L1*L7*MA*cos(Q1)*cos(Q1-Q7)/(L2+Q2)-L7*MB*(2*sin(Q1)*sin(Q1-Q7)+cos(Q1)*(L2+2*Q2)*cos(Q1-Q7)/(L2+Q2))-cos(Q5)*(L5*ME+MF*(L6+2*Q6))*(L7*cos(Q5-Q7)-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/(L6+Q6));
G(2) = -0.5*GRAV*(L1*MA*sin(Q1)*cos(Q1)/(L2+Q2)+L3*MC*sin(Q3)*cos(Q3)/(L4+Q4)+L5*ME*sin(Q5)*cos(Q5)/(L6+Q6)-MB*sin(Q1)*cos(Q1)*(2-(L2+2*Q2)/(L2+Q2))-MD*sin(Q3)*cos(Q3)*(2-(L4+2*Q4)/(L4+Q4))-MF*sin(Q5)*cos(Q5)*(2-(L6+2*Q6)/(L6+Q6)));
G(3) = 0.5*GRAV*(2*MG+2*MB*sin(Q1)^2+2*MD*sin(Q3)^2+2*MF*sin(Q5)^2+cos(Q1)^2*(L1*MA+MB*(L2+2*Q2))/(L2+Q2)+cos(Q3)^2*(L3*MC+MD*(L4+2*Q4))/(L4+Q4)+cos(Q5)^2*(L5*ME+MF*(L6+2*Q6))/(L6+Q6));
GTRAN(1,1) = -L7*sin(Q1-Q7);
GTRAN(1,2) = L7*sin(Q3-Q7);
GTRAN(1,3) = sin(Q5)*(L7*cos(Q7)+L8*sin(Q7)) - L7*sin(Q5-Q7) - cos(Q5)*(L7*sin(Q7)-L8*cos(Q7));
GTRAN(2,1) = cos(Q1);
GTRAN(2,2) = cos(Q3);
GTRAN(2,3) = cos(Q5);
GTRAN(3,1) = sin(Q1);
GTRAN(3,2) = sin(Q3);
GTRAN(3,3) = sin(Q5);



%===========================================================================
function WriteOutput( fileIdentifier, Output )
numberOfOutputQuantities = length( Output );
if numberOfOutputQuantities > 0,
  for i=1:numberOfOutputQuantities,
    fprintf( fileIdentifier, ' %- 14.6E', Output(i) );
  end
  fprintf( fileIdentifier, '\n' );
end



%===========================================================================
% mdlTerminate: Perform end of simulation tasks and set sys=[]
%===========================================================================
function sys = mdlTerminate(T,VAR,u)
FileIdentifier = fopen('all');
fclose( FileIdentifier(1) );
fclose( FileIdentifier(2) );
fprintf( 1, '\n Output is in the files cellsimplifiedmodv4drugged_autolev_2_matlab.i  (i=1,2)\n' );
fprintf( 1, ' The output quantities and associated files are listed in the file cellsimplifiedmodv4drugged_autolev_2_matlab.dir\n' );
fprintf( 1, '\n To load and plot columns 1 and 2 with a solid line and columns 1 and 3 with a dashed line, enter:\n' );
fprintf( 1, '    someName = load( ''cellsimplifiedmodv4drugged_autolev_2_matlab.1'' );\n' );
fprintf( 1, '    plot( someName(:,1), someName(:,2), ''-'', someName(:,1), someName(:,3), ''--'' )\n\n' );
sys = [];



%===========================================================================
% Sfunction: System/Simulink function from standard template
%===========================================================================
function [sys,x0,str,ts] = Sfunction(t,x,u,flag)
switch flag,
  case 0,  [sys,x0,str,ts] = mdlInitializeSizes;    % Initialization of sys, initial state x0, state ordering string str, and sample times ts
  case 1,  sys = mdlDerivatives(t,x,u);             % Calculate the derivatives of continuous states and store them in sys
  case 2,  sys = mdlUpdate(t,x,u);                  % Update discrete states x(n+1) in sys
  case 3,  sys = mdlOutputs(t,x,u);                 % Calculate outputs in sys
  case 4,  sys = mdlGetTimeOfNextVarHit(t,x,u);     % Return next sample time for variable-step in sys
  case 9,  sys = mdlTerminate(t,x,u);               % Perform end of simulation tasks and set sys=[]
  otherwise error(['Unhandled flag = ',num2str(flag)]);
end



%===========================================================================
% mdlInitializeSizes: Return the sizes, initial state VAR, and sample times ts
%===========================================================================
function [sys,VAR,stateOrderingStrings,timeSampling] = mdlInitializeSizes
sizes = simsizes;             % Call simsizes to create a sizes structure
sizes.NumContStates  = 13;    % sys(1) is the number of continuous states
sizes.NumDiscStates  = 0;     % sys(2) is the number of discrete states
sizes.NumOutputs     = 15;    % sys(3) is the number of outputs
sizes.NumInputs      = 0;     % sys(4) is the number of inputs
sizes.DirFeedthrough = 1;     % sys(6) is 1, and allows for the output to be a function of the input
sizes.NumSampleTimes = 1;     % sys(7) is the number of samples times (the number of rows in ts)
sys = simsizes(sizes);        % Convert it to a sizes array
stateOrderingStrings = [];
timeSampling         = [0 0]; % m-by-2 matrix containing the sample times
OpenOutputFilesAndWriteHeadings
VAR = ReadUserInput
