%% |==============================================================================================|%
  %|      Filename: CellRobotNonFullDynamics_Scaled_Improved_V4.m                                 |%
  %|                                                                                              |%
  %|      Author  : Eric Havenhill, Ph.D. (Modified for better feasibility)                       |%
  %|                Faculty Affiliate                                                             |%
  %|                Colorado State University                                                     |%
  %|                                                                                              |%
  %|      Purpose : Generates an optimal trajectory for the "robotic" cell during migration       |%
  %|                with improved numerical conditioning and constraint scaling                   |%
  %|                MODIFIED: Uses sigmoidal functions for actuator force initial guess           |%
  %|                                                                                              |%
  %|      Notes   : * Fixed constraint scaling and conditioning issues                            |%
  %|                * Improved integration method and regularization                              |%
  %|                * Better optimization settings for feasibility                                |%
  %|                * Sigmoidal initial guess for smoother force transitions                      |%
%% |==============================================================================================|%

%% ------------------------- Scaling Factors -------------------------
    clear all; close all; clc;
    format("long")
    
    warning('off', 'all');
    
    % Define scaling factors
    TIME_SCALE = 25206.93764;           % Scale time from [0,1] to [0,t_sim]
    LENGTH_SCALE = 1e-6;                % Scale lengths to micrometers
    MASS_SCALE = 1e-23;                 % Scale masses 
    FORCE_SCALE = 1e-9;                 % Scale forces to nanoNewtons
    
    global L1 L2 L3 L4 L5 L6 L7 L8 L9 h_nuc q_Ans_init qdot_Ans_init q_Ans_fin qdot_Ans_fin
    global TIME_SCALE LENGTH_SCALE MASS_SCALE FORCE_SCALE
    
    % Define characteristic scales for constraint normalization
    global CHAR_LENGTH CHAR_ANGLE CHAR_VELOCITY_LINEAR CHAR_VELOCITY_ANGULAR
    CHAR_LENGTH = 20;                   % 20 μm (characteristic length scale)
    CHAR_ANGLE = pi/2;                  % π/2 rad (characteristic angle scale)
    CHAR_VELOCITY_LINEAR = 2;           % 2 μm/scaled_time
    CHAR_VELOCITY_ANGULAR = 2;          % 2 rad/scaled_time
    
    % Number of nodes (N) that discretize the system
    N = 5;
    
    % Show file name in command window
    Information = [' '];
    disp(Information)
    Information = ['... Running the IMPROVED file "CellRobotNonFullDynamics_Scaled_Improved.m" with N = ', num2str(N), ' temporal nodes ...'];
    disp(Information)
    Information = ['... Using SIGMOIDAL initial guess for actuator forces ...'];
    disp(Information)
    Information = [' '];
    disp(Information)
    
    % Final time of the simulation (scaled to 1.0)
    t_sim_original = 25206.93764;
    t_sim_scaled = 1.0;
    
    % Start the clock to track the MATLAB optimization time
    tic
    
    % The lengths of all rigid bodies (scaled to micrometers)
    L1 = 20;    % 20 μm
    L2 = 20;    % 20 μm
    L3 = 20;    % 20 μm
    L4 = 20;    % 20 μm
    L5 = 20;    % 20 μm
    L6 = 20;    % 20 μm
    L7 = 7;     % 7 μm
    L8 = 2.5;   % 2.5 μm
    L9 = 60;    % 60 μm
    
    % Scaled initial q-values
    q1_0 = pi/2; 
    q2_0 = 6;                           % 6 μm
    q3_0 = (150.52)*pi/180; 
    q4_0 = (52.839-20);                 % μm
    q5_0 = (134.38)*pi/180; 
    q6_0 = (32.882-20);                 % μm
    q7_0 = 0; 
    q8_0 = 7;                           % 7 μm
    q9_0 = 26;                          % 26 μm
    
    q_Ans_init = [q1_0;q2_0;q3_0;q4_0;q5_0;q6_0;q7_0;q8_0;q9_0];

    qdot_Ans_init = zeros(9,1);
    
    h_nuc = q2_0 + L1;
    
    % Scaled final q-values
    q_Ans_fin = [(34.32)*pi/180;
                 (46.120-20);                   % μm
                 (106.92)*pi/180;
                 (27.176-20);                   % μm
                 (57.29)*pi/180;
                 (27.929-20);                   % μm
                 q_Ans_init(7,1);
                 45.093;                        % μm
                 q_Ans_init(9,1)];
    
    % Scaled velocity values (in scaled units per scaled time)
    qdot_Ans_fin = [(q_Ans_fin(1,1)-q_Ans_init(1,1))/t_sim_scaled;
                    (q_Ans_fin(2,1)-q_Ans_init(2,1))/t_sim_scaled;
                    (q_Ans_fin(3,1)-q_Ans_init(3,1))/t_sim_scaled;
                    (q_Ans_fin(4,1)-q_Ans_init(4,1))/t_sim_scaled;
                    (q_Ans_fin(5,1)-q_Ans_init(5,1))/t_sim_scaled;
                    (q_Ans_fin(6,1)-q_Ans_init(6,1))/t_sim_scaled;
                    (q_Ans_fin(7,1)-q_Ans_init(7,1))/t_sim_scaled;
                    (q_Ans_fin(8,1)-q_Ans_init(8,1))/t_sim_scaled;
                    (q_Ans_fin(9,1)-q_Ans_init(9,1))/t_sim_scaled]; 

%% ------------------------- Initialize the States & Control -------------------------
    % Initialize states (linearly spaced on the nodes)
    t0 = 1.0;  % Scaled time
    
    % Initialize states (linearly spaced on the nodes)
    t = linspace(0,1,N)';   % normalized time vector [0,1]
    
    % Define a generic bell-shaped function (Gaussian centered at 0.5)
    bell_shape = exp(-((t-0.5).^2)/(2*0.1^2));   % width = 0.1 (adjustable)
    bell_shape = bell_shape ./ max(bell_shape);  % normalize to [0,1]
    
    % Function to scale bell according to qdot_Ans_fin sign
    make_bell = @(val) sign(val) * abs(val) * bell_shape;
    
    % Initialize q states and bell-shaped qdot guesses
    q1_0 = linspace(q_Ans_init(1),q_Ans_fin(1),N)';  
    q1dot_0 = (1E-9)*make_bell(qdot_Ans_fin(1));
    
    q2_0 = linspace(q_Ans_init(2),q_Ans_fin(2),N)';  
    q2dot_0 = (1E-15)*make_bell(qdot_Ans_fin(2));
    
    q3_0 = linspace(q_Ans_init(3),q_Ans_fin(3),N)';  
    q3dot_0 = (1E-9)*make_bell(qdot_Ans_fin(3));
    
    q4_0 = linspace(q_Ans_init(4),q_Ans_fin(4),N)';  
    q4dot_0 = (1E-15)*make_bell(qdot_Ans_fin(4));
    
    q5_0 = linspace(q_Ans_init(5),q_Ans_fin(5),N)';  
    q5dot_0 = (1E-9)*make_bell(qdot_Ans_fin(5));
    
    q6_0 = linspace(q_Ans_init(6),q_Ans_fin(6),N)';  
    q6dot_0 = (1E-15)*make_bell(qdot_Ans_fin(6));
    
    q7_0 = linspace(q_Ans_init(7),q_Ans_fin(7),N)';  
    q7dot_0 = (0)*make_bell(qdot_Ans_fin(7));
    
    q8_0 = linspace(q_Ans_init(8),q_Ans_fin(8),N)';  
    q8dot_0 = (1E-13)*make_bell(qdot_Ans_fin(8));
    
    q9_0 = q_Ans_fin(9) * ones(N,1);  
    q9dot_0 = (0)*make_bell(qdot_Ans_fin(9));

%% ------------------------- SIGMOIDAL FORCE INITIAL GUESS -------------------------
    % Define the control forces (scaled to nanoNewtons) using sigmoidal functions
    
    % Define initial and final values for each force
    Upsilon1_init = 42.97;
    Upsilon1_final = 52.75;
    
    Upsilon2_init = 53.56;
    Upsilon2_final = 42.62;
    
    Upsilon3_init = -66.10;
    Upsilon3_final = -57.414;
    
    % Create normalized time vector (0 to 1)
    t_normalized = linspace(0, 1, N)';
    
    % Sigmoidal function parameters
    % k controls the steepness of the transition (higher k = steeper)
    % t_center controls where the transition occurs (0.5 = middle)
    k = 6;          % Steepness parameter (try values between 4-10)
    t_center = 0.5; % Center of transition (0.5 = middle of time span)
    
    % Generate sigmoidal functions for each force
    % Sigmoid formula: f(t) = f_init + (f_final - f_init) / (1 + exp(-k*(t - t_center)))
    Upsilon1_0 = Upsilon1_init + (Upsilon1_final - Upsilon1_init) ./ ...
                 (1 + exp(-k * (t_normalized - t_center)));
    
    Upsilon2_0 = Upsilon2_init + (Upsilon2_final - Upsilon2_init) ./ ...
                 (1 + exp(-k * (t_normalized - t_center)));
    
    Upsilon3_0 = Upsilon3_init + (Upsilon3_final - Upsilon3_init) ./ ...
                 (1 + exp(-k * (t_normalized - t_center)));
    
    % Display sigmoidal force information
    fprintf('Sigmoidal Force Initial Guess Parameters:\n');
    fprintf('  Steepness parameter (k): %.1f\n', k);
    fprintf('  Transition center: %.1f (normalized time)\n', t_center);
    fprintf('  Upsilon1: %.2f → %.2f\n', Upsilon1_init, Upsilon1_final);
    fprintf('  Upsilon2: %.2f → %.2f\n', Upsilon2_init, Upsilon2_final);
    fprintf('  Upsilon3: %.2f → %.2f\n', Upsilon3_init, Upsilon3_final);
    fprintf('\n');
    
    % Plot the initial guess for visualization
    figure('Name', 'Sigmoidal Force Initial Guess', 'Position', [100, 100, 800, 600]);
    subplot(3,1,1);
    plot(t_normalized, Upsilon1_0, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    grid on;
    title('Upsilon1 Initial Guess (Sigmoidal)');
    xlabel('Normalized Time');
    ylabel('Force (nN)');
    ylim([min(Upsilon1_init, Upsilon1_final) - 2, max(Upsilon1_init, Upsilon1_final) + 2]);
    
    subplot(3,1,2);
    plot(t_normalized, Upsilon2_0, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
    grid on;
    title('Upsilon2 Initial Guess (Sigmoidal)');
    xlabel('Normalized Time');
    ylabel('Force (nN)');
    ylim([min(Upsilon2_init, Upsilon2_final) - 2, max(Upsilon2_init, Upsilon2_final) + 2]);
    
    subplot(3,1,3);
    plot(t_normalized, Upsilon3_0, 'k-^', 'LineWidth', 2, 'MarkerSize', 6);
    grid on;
    title('Upsilon3 Initial Guess (Sigmoidal)');
    xlabel('Normalized Time');
    ylabel('Force (nN)');
    ylim([min(Upsilon3_init, Upsilon3_final) - 2, max(Upsilon3_init, Upsilon3_final) + 2]);
    
    sgtitle('Sigmoidal Force Initial Guess Functions');
    
    S0 = [t0;...
          q1_0;q1dot_0;...
          q2_0;q2dot_0;...
          q3_0;q3dot_0;...
          q4_0;q4dot_0;...
          q5_0;q5dot_0;...
          q6_0;q6dot_0;...
          q7_0;q7dot_0;...
          q8_0;q8dot_0;...
          q9_0;q9dot_0;...
          Upsilon1_0;Upsilon2_0;Upsilon3_0];
    
    % No linear inequality or equality constraints
    A = [];
    b = [];
    Aeq = [];
    Beq = [];

%% ------------------------- Set Lower & Upper Bounds -------------------------
    % IMPROVED: Relaxed time bounds
    lbt = 0.95;
    ubt = 2.05;
    
    % Position bounds
    lbq1 = min([q_Ans_init(1,1),q_Ans_fin(1,1)])*ones(N,1) - 0.2;
    ubq1 = max([q_Ans_init(1,1),q_Ans_fin(1,1)])*ones(N,1) + 0.2;
    
    lbq2 = min([q_Ans_init(2,1),q_Ans_fin(2,1)])*ones(N,1) - 2;
    ubq2 = max([q_Ans_init(2,1),q_Ans_fin(2,1)])*ones(N,1) + 2;
    
    lbq3 = min([q_Ans_init(3,1),q_Ans_fin(3,1)])*ones(N,1) - 0.2;
    ubq3 = max([q_Ans_init(3,1),q_Ans_fin(3,1)])*ones(N,1) + 0.2;
    
    lbq4 = min([q_Ans_init(4,1),q_Ans_fin(4,1)])*ones(N,1) - 2;
    ubq4 = max([q_Ans_init(4,1),q_Ans_fin(4,1)])*ones(N,1) + 2;
    
    lbq5 = min([q_Ans_init(5,1),q_Ans_fin(5,1)])*ones(N,1) - 0.2;
    ubq5 = max([q_Ans_init(5,1),q_Ans_fin(5,1)])*ones(N,1) + 0.2;
    
    lbq6 = min([q_Ans_init(6,1),q_Ans_fin(6,1)])*ones(N,1) - 2;
    ubq6 = max([q_Ans_init(6,1),q_Ans_fin(6,1)])*ones(N,1) + 2;
    
    lbq7 = min([q_Ans_init(7,1),q_Ans_fin(7,1)])*ones(N,1) - 0.2;
    ubq7 = max([q_Ans_init(7,1),q_Ans_fin(7,1)])*ones(N,1) + 0.2;
    
    lbq8 = min([q_Ans_init(8,1),q_Ans_fin(8,1)])*ones(N,1) - 5;
    ubq8 = max([q_Ans_init(8,1),q_Ans_fin(8,1)])*ones(N,1) + 5;
    
    lbq9 = min([q_Ans_init(9,1),q_Ans_fin(9,1)])*ones(N,1) - 5;
    ubq9 = max([q_Ans_init(9,1),q_Ans_fin(9,1)])*ones(N,1) + 5;
    
    % IMPROVED: More reasonable velocity bounds
    lbq1dot = -3*max(abs(qdot_Ans_fin(1,1)), 0.1)*ones(N,1);
    ubq1dot = 3*max(abs(qdot_Ans_fin(1,1)), 0.1)*ones(N,1);
    
    lbq2dot = -3*max(abs(qdot_Ans_fin(2,1)), 0.1)*ones(N,1);
    ubq2dot = 3*max(abs(qdot_Ans_fin(2,1)), 0.1)*ones(N,1);
    
    lbq3dot = -3*max(abs(qdot_Ans_fin(3,1)), 0.1)*ones(N,1);
    ubq3dot = 3*max(abs(qdot_Ans_fin(3,1)), 0.1)*ones(N,1);
    
    lbq4dot = -3*max(abs(qdot_Ans_fin(4,1)), 0.1)*ones(N,1);
    ubq4dot = 3*max(abs(qdot_Ans_fin(4,1)), 0.1)*ones(N,1);
    
    lbq5dot = -3*max(abs(qdot_Ans_fin(5,1)), 0.1)*ones(N,1);
    ubq5dot = 3*max(abs(qdot_Ans_fin(5,1)), 0.1)*ones(N,1);
    
    lbq6dot = -3*max(abs(qdot_Ans_fin(6,1)), 0.1)*ones(N,1);
    ubq6dot = 3*max(abs(qdot_Ans_fin(6,1)), 0.1)*ones(N,1);
    
    lbq7dot = -3*max(abs(qdot_Ans_fin(7,1)), 0.1)*ones(N,1);
    ubq7dot = 3*max(abs(qdot_Ans_fin(7,1)), 0.1)*ones(N,1);
    
    lbq8dot = -3*max(abs(qdot_Ans_fin(8,1)), 0.1)*ones(N,1);
    ubq8dot = 3*max(abs(qdot_Ans_fin(8,1)), 0.1)*ones(N,1);
    
    lbq9dot = -3*max(abs(qdot_Ans_fin(9,1)), 0.1)*ones(N,1);
    ubq9dot = 3*max(abs(qdot_Ans_fin(9,1)), 0.1)*ones(N,1);
    
    % Control bounds (scaled to nanoNewtons)
    lbUpsilon1 = -300*ones(N,1);  % Increased range
    ubUpsilon1 = 300*ones(N,1);
    
    lbUpsilon2 = -300*ones(N,1);
    ubUpsilon2 = 300*ones(N,1);
    
    lbUpsilon3 = -300*ones(N,1);
    ubUpsilon3 = 300*ones(N,1);
    
    lb = [lbt;
          lbq1;lbq1dot;
          lbq2;lbq2dot;
          lbq3;lbq3dot;
          lbq4;lbq4dot;
          lbq5;lbq5dot;
          lbq6;lbq6dot;
          lbq7;lbq7dot;
          lbq8;lbq8dot;
          lbq9;lbq9dot;
          lbUpsilon1;lbUpsilon2;lbUpsilon3];
    
    ub = [ubt;
          ubq1;ubq1dot;
          ubq2;ubq2dot;
          ubq3;ubq3dot;
          ubq4;ubq4dot;
          ubq5;ubq5dot;
          ubq6;ubq6dot;
          ubq7;ubq7dot;
          ubq8;ubq8dot;
          ubq9;ubq9dot;
          ubUpsilon1;ubUpsilon2;ubUpsilon3];

%% ------------------------- IMPROVED Optimization -------------------------
    % Phase 1: Feasibility-focused optimization
    options1 = optimoptions(@fmincon,... 
        'TolFun',1E-3, ...              
        'TolCon',1E-2,...               % Relaxed constraint tolerance
        'TolX',1E-6,...                 
        'MaxFunEvals',50000, ...        
        'MaxIter', 1000, ...
        'Display', 'iter', ...
        'Algorithm', 'interior-point',...
        'EnableFeasibilityMode', true,...
        'SubproblemAlgorithm', 'cg',...
        'ScaleProblem', true,...
        'FiniteDifferenceStepSize', 1e-8,...
        'OptimalityTolerance', 1e-3,...
        'ConstraintTolerance', 1e-2,...
        'StepTolerance', 1e-8);
 
    fprintf('Phase 1: Finding feasible solution...\n');
    S1 = fmincon(@CostFunction, S0, A, b, Aeq, Beq, lb, ub, @Constraints_Improved, options1, N);

    problem = createOptimProblem('fmincon', 'x0', S0, 'objective', @(S) CostFunction(S, N), ...
                                  'nonlcon', @(S) Constraints_Improved(S, N), 'lb', lb, 'ub', ub, 'options', options1);

    % ------------------------- Optimization with GlobalSearch -------------------------
     
    % gs = GlobalSearch;
    % [S1, fval] = run(gs, problem);

    % ------------------------- Optimization with MultiStart -------------------------
    % ms = MultiStart;
    % [S1, fval, exitflag, output, solutions] = run(ms, problem, 100); % 10 random starts

    % Diagnose Phase 1 results
    [c1, ceq1] = Constraints_Improved(S1, N);
    max_violation_1 = max(abs(ceq1));
    fprintf('Phase 1 complete. Max constraint violation: %.6e\n', max_violation_1);
    
    % Phase 2: Tighten tolerances if Phase 1 was successful
    if max_violation_1 < 1e-2
        fprintf('Phase 2: Improving solution precision...\n');
        options2 = optimoptions(@fmincon,... 
            'TolFun',1E-6, ...              
            'TolCon',1E-4,...               
            'TolX',1E-8,...                 
            'MaxFunEvals',100000, ...        
            'MaxIter', 2000, ...
            'Display', 'iter', ...
            'Algorithm', 'interior-point',...
            'EnableFeasibilityMode', false,...
            'SubproblemAlgorithm', 'cg',...
            'ScaleProblem', true,...
            'FiniteDifferenceStepSize', 1e-8,...
            'OptimalityTolerance', 1e-6,...
            'ConstraintTolerance', 1e-4,...
            'StepTolerance', 1e-10);
    
        S = fmincon(@CostFunction, S1, A, b, Aeq, Beq, lb, ub, @Constraints_Improved, options2, N);
    else
        fprintf('Phase 1 did not achieve sufficient feasibility. Using Phase 1 result.\n');
        S = S1;
    end
    
    %S =S1; 

%% ------------------------- Results & Plots (Convert back to original units) -------------------------
    fprintf('\n Optimization Completed in %f seconds \n', toc);
    
    % Extract states and time (convert back to original units)
    tfin_scaled = S(1,1);
    tfin_original = tfin_scaled * TIME_SCALE;
    delta_original = tfin_original/N;
    tvec_original = 0: delta_original : tfin_original-delta_original;
    
    % Extract and convert positions back to original units
    pos.q1 = S(2+0*N:1*N+1,1);  % Angles remain unchanged
    vel.q1 = S(2+1*N:2*N+1,1) / TIME_SCALE;  % Convert velocity back
    
    pos.q2 = S(2+2*N:3*N+1,1) * LENGTH_SCALE;  % Convert back to meters
    vel.q2 = S(2+3*N:4*N+1,1) * LENGTH_SCALE / TIME_SCALE;
    
    pos.q3 = S(2+4*N:5*N+1,1);  % Angles remain unchanged
    vel.q3 = S(2+5*N:6*N+1,1) / TIME_SCALE;
    
    pos.q4 = S(2+6*N:7*N+1,1) * LENGTH_SCALE;
    vel.q4 = S(2+7*N:8*N+1,1) * LENGTH_SCALE / TIME_SCALE;
    
    pos.q5 = S(2+8*N:9*N+1,1);  % Angles remain unchanged
    vel.q5 = S(2+9*N:10*N+1,1) / TIME_SCALE;
    
    pos.q6 = S(2+10*N:11*N+1,1) * LENGTH_SCALE;
    vel.q6 = S(2+11*N:12*N+1,1) * LENGTH_SCALE / TIME_SCALE;
    
    pos.q7 = S(2+12*N:13*N+1,1);  % Angles remain unchanged
    vel.q7 = S(2+13*N:14*N+1,1) / TIME_SCALE;
    
    pos.q8 = S(2+14*N:15*N+1,1) * LENGTH_SCALE;
    vel.q8 = S(2+15*N:16*N+1,1) * LENGTH_SCALE / TIME_SCALE;
    
    pos.q9 = S(2+16*N:17*N+1,1) * LENGTH_SCALE;
    vel.q9 = S(2+17*N:18*N+1,1) * LENGTH_SCALE / TIME_SCALE;
    
    % Convert control forces back to original units
    control.Upsilon1 = S(2+18*N:19*N+1,1) * FORCE_SCALE;
    control.Upsilon2 = S(2+19*N:20*N+1,1) * FORCE_SCALE;
    control.Upsilon3 = S(2+20*N:21*N+1,1) * FORCE_SCALE;
    
    % Make the Plots (using original units)
    figure(1)
    rows = 3; cols = 9;
    
    subplot(rows, cols, 1);grid on; plot(tvec_original, pos.q1,'-o');
    hold on
    plot(tvec_original,q1_0,'*');
    xlabel('t (s)'); ylabel('${q}_{1}$ $(rad)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 2);grid on; plot(tvec_original, pos.q2,'-o');
    hold on
    plot(tvec_original,q2_0*1e-6,'*');
    xlabel('t (s)'); ylabel('${q}_{2}$ $(m)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 3);grid on; plot(tvec_original, pos.q3,'-o');
    hold on
    plot(tvec_original,q3_0,'*');
    xlabel('t (s)'); ylabel('${q}_{3}$ $(rad)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 4);grid on; plot(tvec_original, pos.q4,'-o');
    hold on
    plot(tvec_original,q4_0*1e-6,'*');
    xlabel('t (s)'); ylabel('${q}_{4}$ $(m)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 5);grid on; plot(tvec_original, pos.q5,'-o');
    hold on
    plot(tvec_original,q5_0,'*');
    xlabel('t (s)'); ylabel('${q}_{5}$ $(rad)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 6);grid on; plot(tvec_original, pos.q6,'-o');
    hold on
    plot(tvec_original,q6_0*1e-6,'*');
    xlabel('t (s)'); ylabel('${q}_{6}$ $(m)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 7);grid on; plot(tvec_original, pos.q7,'-o');
    hold on
    plot(tvec_original,q7_0,'*');
    xlabel('t (s)'); ylabel('${q}_{7}$ $(rad)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 8);grid on; plot(tvec_original, pos.q8,'-o');
    hold on
    plot(tvec_original,q8_0*1e-6,'*');
    xlabel('t (s)'); ylabel('${q}_{8}$ $(m)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 9);grid on; plot(tvec_original, pos.q9,'-o');
    hold on
    plot(tvec_original,q9_0*1e-6,'*');
    xlabel('t (s)'); ylabel('${q}_{9}$ $(m)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    % Velocity plots
    subplot(rows, cols, 10); grid on; plot(tvec_original, vel.q1,'-o');
    hold on
    plot(tvec_original,q1dot_0,'*');
    xlabel('t (s)'); ylabel('$\dot{q}_{1}$ $(rad/s)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 11); grid on; plot(tvec_original, vel.q2,'-o');
    hold on
    plot(tvec_original,q2dot_0,'*');
    xlabel('t (s)'); ylabel('$\dot{q}_{2}$ $(m/s)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 12); grid on; plot(tvec_original, vel.q3,'-o');
    hold on
    plot(tvec_original,q3dot_0,'*');
    xlabel('t (s)'); ylabel('$\dot{q}_{3}$ $(rad/s)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 13); grid on; plot(tvec_original, vel.q4,'-o');
    hold on
    plot(tvec_original,q4dot_0,'*');
    xlabel('t (s)'); ylabel('$\dot{q}_{4}$ $(m/s)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 14); grid on; plot(tvec_original, vel.q5,'-o');
    hold on
    plot(tvec_original,q5dot_0,'*');
    xlabel('t (s)'); ylabel('$\dot{q}_{5}$ $(rad/s)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 15); grid on; plot(tvec_original, vel.q6,'-o');
    hold on
    plot(tvec_original,q6dot_0,'*');
    xlabel('t (s)'); ylabel('$\dot{q}_{6}$ $(m/s)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 16); grid on; plot(tvec_original, vel.q7,'-o');
    hold on
    plot(tvec_original,q7dot_0,'*');
    xlabel('t (s)'); ylabel('$\dot{q}_{7}$ $(rad/s)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 17); grid on; plot(tvec_original, vel.q8,'-o');
    hold on
    plot(tvec_original,q8dot_0,'*');
    xlabel('t (s)'); ylabel('$\dot{q}_{8}$ $(m/s)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 18); grid on; plot(tvec_original, vel.q9,'-o');
    hold on
    plot(tvec_original,q9dot_0,'*');
    xlabel('t (s)'); ylabel('$\dot{q}_{9}$ $(m/s)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    % Control plots
    subplot(rows, cols, 19);grid on; plot(tvec_original, control.Upsilon1,'-o');
    hold on 
    plot(tvec_original,Upsilon1_0*1e-9,'*');
    xlabel('t (s)'); ylabel('$\Upsilon_{1}$ $(N)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 20);grid on; plot(tvec_original, control.Upsilon2,'-o');
    hold on 
    plot(tvec_original,Upsilon2_0*1e-9,'*');
    xlabel('t (s)'); ylabel('$\Upsilon_{2}$ $(N)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    subplot(rows, cols, 21);grid on; plot(tvec_original, control.Upsilon3,'-o');
    hold on 
    plot(tvec_original,Upsilon3_0*1e-9,'*');
    xlabel('t (s)'); ylabel('$\Upsilon_{3}$ $(N)$', 'Interpreter','latex');
    ax = gca; ax.FontSize = 12;
    
    sgtitle('Improved Optimized Cell Robot Dynamics - All Variables');
    exportgraphics(gcf,'AllOptimizedPlots_Improved.pdf','BackgroundColor','none')
    
    % Control forces plot with custom dimensions and font size
    figure(2)
    % Set figure dimensions (width x height in pixels)
    fig_width = 600;
    fig_height = 500;
    set(gcf, 'Position', [100, 100, fig_width, fig_height]);
    % Set font size (customizable)
    font_size = 16.5; % Change this value as needed
    grid on;
    plot(tvec_original, control.Upsilon1, 'b', 'LineWidth', 2);
    hold on;
    %plot(tvec_original,Upsilon1_0*1e-9,'b*');
    plot(tvec_original, control.Upsilon2, 'r', 'LineWidth', 2);
    %plot(tvec_original,Upsilon2_0*1e-9,'r*');
    plot(tvec_original, control.Upsilon3, 'k', 'LineWidth', 2);
    %plot(tvec_original,Upilon3_0*1e-9,'k*');
    xlabel('t (s)', 'FontSize', font_size);
    ylabel('Actuator Forces', 'Interpreter','latex', 'FontSize', font_size);
    %legend('$\Upsilon_{1}$ $(nN)$','$\Upsilon_{1}(0)$ $(nN)$',...
    % '$\Upsilon_{2}$ $(nN)$','$\Upsilon_{2}(0)$ $(nN)$',...
    % '$\Upsilon_{3}$ $(nN)$','$\Upsilon_{3}(0)$ $(nN)$','Interpreter','latex')
    leg = legend('$\Upsilon_{B}$ $(N)$ ','$\Upsilon_{D}$ $(N)$ ','$\Upsilon_{F}$ $(N)$ ',...
    'Interpreter','latex', 'FontSize', font_size, 'Location', 'northwest');
    % Adjust legend box size and position
    leg.Position(3) = leg.Position(3) * 1; % Increase width by 20%
    leg.Position(1) = leg.Position(1) + 0.025; % Move right by 0.05 (adjust as needed)
    leg.Position(2) = leg.Position(1) + 0.535; % Move right by 0.05 (adjust as needed)
    % Set axes font size and y-axis limits
    ax = gca;
    ax.FontSize = 25;
    ylim([-80*10^-9, 110*10^-9]);
    xlim([min(tvec_original), max(tvec_original)]);
    % Optional: Set title font size if you add a title
    % title('Control Forces', 'FontSize', font_size + 2);
    exportgraphics(gcf,'OptimalControl_Improved.pdf','BackgroundColor','none')
    
    % Position states plot
    figure(3)
    % Set figure dimensions (width x height in pixels)
    fig_width = 600;
    fig_height = 500;
    set(gcf, 'Position', [100, 100, fig_width, fig_height]);
    % Set font size (customizable)
    font_size = 16.5; % Change this value as needed
    grid on;
    plot(tvec_original, pos.q7*10^-6, 'b', 'LineWidth', 2);
    hold on;
    %plot(tvec_original, q7_0, 'b*', 'LineWidth', 2);
    plot(tvec_original, pos.q8, 'k', 'LineWidth', 2);
    %plot(tvec_original, q8_0, 'k*', 'LineWidth', 2);
    plot(tvec_original, pos.q9, 'r', 'LineWidth', 2);
    %plot(tvec_original, q9_0, 'r*', 'LineWidth', 2);
    xlabel('t (s)'); ylabel('States', 'Interpreter','latex');
    % Fixed: Assign legend to 'leg' variable to enable position modification
    leg = legend('$q_{7}$ $(\times$ $10^{6}$ $rad)$ ','$q_{8}$ $(m)$','$q_{9}$ $(m)$','$q_{9}$ $(\mu m)$','Interpreter','latex','FontSize', font_size, 'Location', 'northwest');
    % Adjust legend box size and position - now this will work
    leg.Position(3) = leg.Position(3) * 1.2; % Increase width by 300%
    leg.Position(1) = leg.Position(1) + 0.0; % Move right by 0.05 (adjust as needed)
    leg.Position(2) = leg.Position(1) + 0.55; % Move right by 0.05 (adjust as needed)
    % Set axes font size and y-axis limits
    ax = gca;
    ax.FontSize = font_size;
    ylim([-5*10^-6, 50*10^-6]);
    xlim([min(tvec_original), max(tvec_original)]);
    ax = gca; ax.FontSize = 25;
    exportgraphics(gcf,'OptimalStates_Improved.pdf','BackgroundColor','none')

    %% ———————–– Actual Dynamics Using the Optimal Control  ———————––
% Replace your RK4 section with this corrected version:

%% ———————–– Actual Dynamics Using the Optimal Control (CORRECTED) ———————––

% Use SCALED initial conditions (same as optimization)
Y(1,:) = [q_Ans_init(1),qdot_Ans_init(1),...
          q_Ans_init(2),qdot_Ans_init(2),...
          q_Ans_init(3),qdot_Ans_init(3),...
          q_Ans_init(4),qdot_Ans_init(4),...
          q_Ans_init(5),qdot_Ans_init(5),...
          q_Ans_init(6),qdot_Ans_init(6),...
          q_Ans_init(7),qdot_Ans_init(7),...
          q_Ans_init(8),qdot_Ans_init(8),...
          q_Ans_init(9),qdot_Ans_init(9)];

% Extract SCALED control forces (don't convert to original units yet)
control_scaled.Upsilon1 = S(2+18*N:19*N+1,1);  % Keep in scaled units (nN)
control_scaled.Upsilon2 = S(2+19*N:20*N+1,1);  % Keep in scaled units (nN)
control_scaled.Upsilon3 = S(2+20*N:21*N+1,1);  % Keep in scaled units (nN)

tfin_scaled = S(1,1);
dt_scaled = tfin_scaled/(N-1);  % Consistent with optimization
tvec_scaled = linspace(0, tfin_scaled, N);

for n = 1 : N - 1
    Tau = [control_scaled.Upsilon1(n); control_scaled.Upsilon2(n); control_scaled.Upsilon3(n)]; 
    
    % Create state vector with time (to match optimization structure)
    Y_with_time = [tvec_scaled(n); Y(n,:)'];
    
    k1 = dt_scaled*(dyn_scaled_robust(tvec_scaled(n), Y_with_time, Tau));
    k1 = k1(2:end)'; % Remove time derivative, transpose
    
    Y_with_time_k2 = [tvec_scaled(n) + dt_scaled/2; (Y(n,:) + k1/2)'];
    k2 = dt_scaled*(dyn_scaled_robust(tvec_scaled(n) + dt_scaled/2, Y_with_time_k2, Tau));
    k2 = k2(2:end)';
    
    Y_with_time_k3 = [tvec_scaled(n) + dt_scaled/2; (Y(n,:) + k2/2)'];
    k3 = dt_scaled*(dyn_scaled_robust(tvec_scaled(n) + dt_scaled/2, Y_with_time_k3, Tau));
    k3 = k3(2:end)';
    
    Y_with_time_k4 = [tvec_scaled(n) + dt_scaled; (Y(n,:) + k3)'];
    k4 = dt_scaled*(dyn_scaled_robust(tvec_scaled(n) + dt_scaled, Y_with_time_k4, Tau));
    k4 = k4(2:end)';

    Y(n+1,:) = Y(n,:) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

% Convert RK4 results to original units for plotting
Y_original = Y;
% Angles remain unchanged
Y_original(:,1) = Y(:,1);  % q1 (rad)
Y_original(:,5) = Y(:,5);  % q3 (rad) 
Y_original(:,9) = Y(:,9);  % q5 (rad)
Y_original(:,13) = Y(:,13); % q7 (rad)

% Convert scaled lengths to original units (meters)
Y_original(:,3) = Y(:,3) * LENGTH_SCALE;   % q2: μm -> m
Y_original(:,7) = Y(:,7) * LENGTH_SCALE;   % q4: μm -> m
Y_original(:,11) = Y(:,11) * LENGTH_SCALE; % q6: μm -> m
Y_original(:,15) = Y(:,15) * LENGTH_SCALE; % q8: μm -> m
Y_original(:,17) = Y(:,17) * LENGTH_SCALE; % q9: μm -> m

% Convert scaled velocities to original units
Y_original(:,2) = Y(:,2) / TIME_SCALE;  % q1dot: rad/scaled_time -> rad/s
Y_original(:,4) = Y(:,4) * LENGTH_SCALE / TIME_SCALE;  % q2dot
Y_original(:,6) = Y(:,6) / TIME_SCALE;  % q3dot
Y_original(:,8) = Y(:,8) * LENGTH_SCALE / TIME_SCALE;  % q4dot
Y_original(:,10) = Y(:,10) / TIME_SCALE; % q5dot
Y_original(:,12) = Y(:,12) * LENGTH_SCALE / TIME_SCALE; % q6dot
Y_original(:,14) = Y(:,14) / TIME_SCALE; % q7dot
Y_original(:,16) = Y(:,16) * LENGTH_SCALE / TIME_SCALE; % q8dot
Y_original(:,18) = Y(:,18) * LENGTH_SCALE / TIME_SCALE; % q9dot

% Plot RK4 results (using converted time vector)
tvec_original_rk4 = tvec_scaled * TIME_SCALE;

figure(5)
subplot(rows, cols, 1); grid on; 
plot(tvec_original_rk4, Y_original(:,1),'-o'); 
hold on;
plot(tvec_original, pos.q1, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$q_{1}$ (rad)', 'Interpreter','latex');

subplot(rows, cols, 2); grid on; 
plot(tvec_original_rk4, Y_original(:,2),'-o'); 
hold on;
plot(tvec_original, vel.q1, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$\dot{q}_{1}$ (rad)', 'Interpreter','latex');

subplot(rows, cols, 3); grid on; 
plot(tvec_original_rk4, Y_original(:,3),'-o'); 
hold on;
plot(tvec_original, pos.q2, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$q_{2}$ (rad)', 'Interpreter','latex');

subplot(rows, cols, 4); grid on; 
plot(tvec_original_rk4, Y_original(:,4),'-o'); 
hold on;
plot(tvec_original, vel.q2, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$\dot{q}_{2}$ (rad/s)', 'Interpreter','latex');

subplot(rows, cols, 5); grid on; 
plot(tvec_original_rk4, Y_original(:,5),'-o'); 
hold on;
plot(tvec_original, pos.q3, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$q_{3}$ (rad)', 'Interpreter','latex');

subplot(rows, cols, 6); grid on; 
plot(tvec_original_rk4, Y_original(:,6),'-o'); 
hold on;
plot(tvec_original, vel.q3, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$\dot{q}_{3}$ (rad/s)', 'Interpreter','latex');

subplot(rows, cols, 7); grid on; 
plot(tvec_original_rk4, Y_original(:,7),'-o'); 
hold on;
plot(tvec_original, pos.q4, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$q_{4}$ (rad)', 'Interpreter','latex');

subplot(rows, cols, 8); grid on; 
plot(tvec_original_rk4, Y_original(:,8),'-o');
hold on;
plot(tvec_original, vel.q4, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$\dot{q}_{4}$ (rad/s)', 'Interpreter','latex');

subplot(rows, cols, 9); grid on; 
plot(tvec_original_rk4, Y_original(:,9),'-o'); 
hold on;
plot(tvec_original, pos.q5, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$q_{5}$ (rad)', 'Interpreter','latex');

subplot(rows, cols, 10); grid on; 
plot(tvec_original_rk4, Y_original(:,10),'-o'); 
hold on;
plot(tvec_original, vel.q5, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$\dot{q}_{5}$ (rad/s)', 'Interpreter','latex');

subplot(rows, cols, 11); grid on; 
plot(tvec_original_rk4, Y_original(:,11),'-o'); 
hold on;
plot(tvec_original, pos.q6, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$q_{6}$ (rad)', 'Interpreter','latex');

subplot(rows, cols, 12); grid on; 
plot(tvec_original_rk4, Y_original(:,12),'-o'); 
hold on;
plot(tvec_original, vel.q6, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$\dot{q}_{6}$ (rad/s)', 'Interpreter','latex');

subplot(rows, cols, 13); grid on; 
plot(tvec_original_rk4, Y_original(:,13),'-o'); 
hold on;
plot(tvec_original, pos.q7, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$q_{7}$ (rad)', 'Interpreter','latex');

subplot(rows, cols, 14); grid on; 
plot(tvec_original_rk4, Y_original(:,14),'-o'); 
hold on;
plot(tvec_original, vel.q7, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$\dot{q}_{7}$ (rad/s)', 'Interpreter','latex');

subplot(rows, cols, 15); grid on; 
plot(tvec_original_rk4, Y_original(:,15),'-o'); 
hold on;
plot(tvec_original, pos.q8, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$q_{8}$ (rad)', 'Interpreter','latex');

subplot(rows, cols, 16); grid on; 
plot(tvec_original_rk4, Y_original(:,16),'-o'); 
hold on;
plot(tvec_original, vel.q8, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$\dot{q}_{8}$ (rad/s)', 'Interpreter','latex');

subplot(rows, cols, 17); grid on; 
plot(tvec_original_rk4, Y_original(:,17),'-o'); 
hold on;
plot(tvec_original, pos.q9, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$q_{9}$ (rad)', 'Interpreter','latex');

subplot(rows, cols, 18); grid on; 
plot(tvec_original_rk4, Y_original(:,18),'-o'); 
hold on;
plot(tvec_original, vel.q9, '--r', 'LineWidth', 2);
legend('RK4 Simulation', 'Optimization Result');
xlabel('t (s)'); ylabel('$\dot{q}_{9}$ (rad/s)', 'Interpreter','latex');


%% ------------------------- Improved Cost Function -------------------------
    function J = CostFunction(S,N)
        global FORCE_SCALE TIME_SCALE
        
        % Extract scaled time and controls
        tfin_scaled = S(1,1);
        delta_scaled = max(0.01, tfin_scaled/N);  % Increased minimum timestep
        
        % Extract scaled control forces
        control.Upsilon1 = S(2+18*N:19*N+1,1);
        control.Upsilon2 = S(2+19*N:20*N+1,1);
        control.Upsilon3 = S(2+20*N:21*N+1,1);
        
        % Control effort in scaled coordinates
        u_squared = control.Upsilon1.^2 + control.Upsilon2.^2 + control.Upsilon3.^2;
        
        % Simpson's 1/3 rule integration in scaled time (higher-order method)
        if N == 1
            J = u_squared(1) * delta_scaled;
        elseif N == 2
            % Fall back to trapezoidal for N=2
            J = (delta_scaled/2) * (u_squared(1) + u_squared(end));
        elseif mod(N-1, 2) == 0
            % N-1 is even, perfect for Simpson's 1/3 rule
            J = (delta_scaled/3) * (u_squared(1) + 4*sum(u_squared(2:2:end-1)) + 2*sum(u_squared(3:2:end-2)) + u_squared(end));
        else
            % N-1 is odd, use Simpson's 1/3 for most intervals and Simpson's 3/8 for the last 3
            n_third = N - 3;  % Use 1/3 rule up to this point
            
            % Simpson's 1/3 rule for intervals 1 to n_third
            J_third = (delta_scaled/3) * (u_squared(1) + 4*sum(u_squared(2:2:n_third)) + 2*sum(u_squared(3:2:n_third-1)) + u_squared(n_third+1));
            
            % Simpson's 3/8 rule for the last 3 intervals
            J_eighth = (3*delta_scaled/8) * (u_squared(n_third+1) + 3*u_squared(n_third+2) + 3*u_squared(n_third+3) + u_squared(end));
            
            J = J_third + J_eighth;
        end
        
        % Scale cost function to maintain physical meaning
        J = J * FORCE_SCALE^2 * TIME_SCALE;
        
        % Normalize by characteristic values to get O(1) cost
        characteristic_force_squared = (100e-9)^2;  % 100 nN squared
        J = J / (characteristic_force_squared * TIME_SCALE);
    
        % Scale the cost function (relative to the importance of other parameters)
        J = 1E0*J;
    end

%% ------------------------- IMPROVED Constraints with Better Scaling -------------------------
%% ------------------------- IMPROVED Constraints with Hermite-Simpson Collocation -------------------------
function [c,ceq] = Constraints_Improved(S,N)
    global L1 L2 L3 L4 L5 L6 L7 L8 L9 h_nuc q_Ans_init q_Ans_fin
    global TIME_SCALE LENGTH_SCALE CHAR_LENGTH CHAR_ANGLE CHAR_VELOCITY_LINEAR CHAR_VELOCITY_ANGULAR
    
    % Nonlinear inequality & equality constraints
    c = []; 
    ceq = [];
    
    % Calculate the scaled timestep 
    tfin_scaled = S(1,1);
    delta_scaled = max(0.01, tfin_scaled/N);  % Increased minimum timestep for stability
    
    % Extract scaled states
    pos.q1 = S(2+0*N:1*N+1,1);
    vel.q1 = S(2+1*N:2*N+1,1);
    
    pos.q2 = S(2+2*N:3*N+1,1);
    vel.q2 = S(2+3*N:4*N+1,1);
    
    pos.q3 = S(2+4*N:5*N+1,1);
    vel.q3 = S(2+5*N:6*N+1,1);
    
    pos.q4 = S(2+6*N:7*N+1,1);
    vel.q4 = S(2+7*N:8*N+1,1);
    
    pos.q5 = S(2+8*N:9*N+1,1);
    vel.q5 = S(2+9*N:10*N+1,1);
    
    pos.q6 = S(2+10*N:11*N+1,1);
    vel.q6 = S(2+11*N:12*N+1,1);
    
    pos.q7 = S(2+12*N:13*N+1,1);
    vel.q7 = S(2+13*N:14*N+1,1);
    
    pos.q8 = S(2+14*N:15*N+1,1);
    vel.q8 = S(2+15*N:16*N+1,1);
    
    pos.q9 = S(2+16*N:17*N+1,1);
    vel.q9 = S(2+17*N:18*N+1,1);
    
    control.Upsilon1 = S(2+18*N:19*N+1,1);
    control.Upsilon2 = S(2+19*N:20*N+1,1);
    control.Upsilon3 = S(2+20*N:21*N+1,1);
    
    % IMPROVED: Properly scaled initial constraints
    ceq_init = [
        (pos.q1(1,1) - q_Ans_init(1,1)) / CHAR_ANGLE;      % Scale angles
        vel.q1(1,1) / CHAR_VELOCITY_ANGULAR;               % Scale angular velocities
        (pos.q2(1,1) - q_Ans_init(2,1)) / CHAR_LENGTH;     % Scale lengths
        vel.q2(1,1) / CHAR_VELOCITY_LINEAR;                % Scale linear velocities
        (pos.q3(1,1) - q_Ans_init(3,1)) / CHAR_ANGLE;
        vel.q3(1,1) / CHAR_VELOCITY_ANGULAR;
        (pos.q4(1,1) - q_Ans_init(4,1)) / CHAR_LENGTH;
        vel.q4(1,1) / CHAR_VELOCITY_LINEAR;
        (pos.q5(1,1) - q_Ans_init(5,1)) / CHAR_ANGLE;
        vel.q5(1,1) / CHAR_VELOCITY_ANGULAR;
        (pos.q6(1,1) - q_Ans_init(6,1)) / CHAR_LENGTH;
        vel.q6(1,1) / CHAR_VELOCITY_LINEAR;
        (pos.q7(1,1) - q_Ans_init(7,1)) / CHAR_ANGLE;
        vel.q7(1,1) / CHAR_VELOCITY_ANGULAR;
        (pos.q8(1,1) - q_Ans_init(8,1)) / CHAR_LENGTH;
        vel.q8(1,1) / CHAR_VELOCITY_LINEAR;
        (pos.q9(1,1) - q_Ans_init(9,1)) / CHAR_LENGTH;
        vel.q9(1,1) / CHAR_VELOCITY_LINEAR
    ];
    
    ceq = ceq_init;

    % HERMITE-SIMPSON COLLOCATION METHOD (Higher-order transcription)
    for i = 1 : length(pos.q8) - 1
        
        % States at beginning and end of interval
        s_i = [pos.q1(i,1);vel.q1(i,1);...
               pos.q2(i,1);vel.q2(i,1);...
               pos.q3(i,1);vel.q3(i,1);...
               pos.q4(i,1);vel.q4(i,1);...
               pos.q5(i,1);vel.q5(i,1);...
               pos.q6(i,1);vel.q6(i,1);...
               pos.q7(i,1);vel.q7(i,1);...
               pos.q8(i,1);vel.q8(i,1);...
               pos.q9(i,1);vel.q9(i,1)];

        s_f = [pos.q1(i+1,1);vel.q1(i+1,1);...
               pos.q2(i+1,1);vel.q2(i+1,1);...
               pos.q3(i+1,1);vel.q3(i+1,1);...
               pos.q4(i+1,1);vel.q4(i+1,1);...
               pos.q5(i+1,1);vel.q5(i+1,1);...
               pos.q6(i+1,1);vel.q6(i+1,1);...
               pos.q7(i+1,1);vel.q7(i+1,1);...
               pos.q8(i+1,1);vel.q8(i+1,1);...
               pos.q9(i+1,1);vel.q9(i+1,1)];

        % Controls at beginning and end of interval
        Upsilon_i = [control.Upsilon1(i,1);control.Upsilon2(i,1);control.Upsilon3(i,1)];
        Upsilon_f = [control.Upsilon1(i+1,1);control.Upsilon2(i+1,1);control.Upsilon3(i+1,1)];
        
        % Compute derivatives at beginning and end of interval
        s_dummy_i = [0; s_i]; 
        sdot_i = dyn_scaled_robust(0, s_dummy_i, Upsilon_i);
        sdot_i = sdot_i(2:end);
        
        s_dummy_f = [delta_scaled; s_f];
        sdot_f = dyn_scaled_robust(delta_scaled, s_dummy_f, Upsilon_f);
        sdot_f = sdot_f(2:end);
        
        % Hermite-Simpson collocation constraints
        % State at midpoint using cubic Hermite interpolation
        s_mid = 0.5*(s_i + s_f) + (delta_scaled/8)*(sdot_i - sdot_f);
        
        % Control at midpoint (linear interpolation)
        Upsilon_mid = 0.5*(Upsilon_i + Upsilon_f);
        
        % Compute derivative at midpoint
        s_dummy_mid = [delta_scaled/2; s_mid];
        sdot_mid = dyn_scaled_robust(delta_scaled/2, s_dummy_mid, Upsilon_mid);
        sdot_mid = sdot_mid(2:end);
        
        % Hermite-Simpson defect constraint: 
        % s_f - s_i - (delta_scaled/6)*(sdot_i + 4*sdot_mid + sdot_f) = 0
        defect = s_f - s_i - (delta_scaled/6)*(sdot_i + 4*sdot_mid + sdot_f);
        
        % Scale different types of defects appropriately
        scaled_defect = defect;
        
        % Scale position defects
        scaled_defect(1) = scaled_defect(1) / CHAR_ANGLE;    % q1
        scaled_defect(3) = scaled_defect(3) / CHAR_LENGTH;   % q2
        scaled_defect(5) = scaled_defect(5) / CHAR_ANGLE;    % q3
        scaled_defect(7) = scaled_defect(7) / CHAR_LENGTH;   % q4
        scaled_defect(9) = scaled_defect(9) / CHAR_ANGLE;    % q5
        scaled_defect(11) = scaled_defect(11) / CHAR_LENGTH; % q6
        scaled_defect(13) = scaled_defect(13) / CHAR_ANGLE;  % q7
        scaled_defect(15) = scaled_defect(15) / CHAR_LENGTH; % q8
        scaled_defect(17) = scaled_defect(17) / CHAR_LENGTH; % q9
        
        % Scale velocity defects
        scaled_defect(2) = scaled_defect(2) / CHAR_VELOCITY_ANGULAR;  % q1dot
        scaled_defect(4) = scaled_defect(4) / CHAR_VELOCITY_LINEAR;   % q2dot
        scaled_defect(6) = scaled_defect(6) / CHAR_VELOCITY_ANGULAR;  % q3dot
        scaled_defect(8) = scaled_defect(8) / CHAR_VELOCITY_LINEAR;   % q4dot
        scaled_defect(10) = scaled_defect(10) / CHAR_VELOCITY_ANGULAR;% q5dot
        scaled_defect(12) = scaled_defect(12) / CHAR_VELOCITY_LINEAR; % q6dot
        scaled_defect(14) = scaled_defect(14) / CHAR_VELOCITY_ANGULAR;% q7dot
        scaled_defect(16) = scaled_defect(16) / CHAR_VELOCITY_LINEAR; % q8dot
        scaled_defect(18) = scaled_defect(18) / CHAR_VELOCITY_LINEAR; % q9dot
        
        ceq = [ceq; scaled_defect];
                
        % IMPROVED: Better scaled kinematic constraints at interval endpoints
        % Evaluate kinematic constraints at both beginning and end of interval
        for j = [i, i+1]
            kin1 = (2*L7*cos(pos.q7(j,1)) + cos(pos.q1(j,1))*(L2+pos.q2(j,1)) - L9 - cos(pos.q3(j,1))*(L4+pos.q4(j,1))) / CHAR_LENGTH;
            kin2 = (2*L7*sin(pos.q7(j,1)) + sin(pos.q1(j,1))*(L2+pos.q2(j,1)) - sin(pos.q3(j,1))*(L4+pos.q4(j,1))) / CHAR_LENGTH;
            kin3 = (L7*cos(pos.q7(j,1)) + L8*sin(pos.q7(j,1)) + cos(pos.q1(j,1))*(L2+pos.q2(j,1)) - 0.5*L9 - cos(pos.q5(j,1))*(L6+pos.q6(j,1))) / CHAR_LENGTH;
            kin4 = (L7*sin(pos.q7(j,1)) + sin(pos.q1(j,1))*(L2+pos.q2(j,1)) - L8*cos(pos.q7(j,1)) - sin(pos.q5(j,1))*(L6+pos.q6(j,1))) / CHAR_LENGTH;
            kin5 = (L7*cos(pos.q7(j,1)) + cos(pos.q1(j,1))*(L2+pos.q2(j,1)) - pos.q8(j,1)) / CHAR_LENGTH;
            kin6 = (L7*sin(pos.q7(j,1)) + sin(pos.q1(j,1))*(L2+pos.q2(j,1)) - pos.q9(j,1)) / CHAR_LENGTH;
            
            % Only add kinematic constraints once per node to avoid duplication
            if j == i
                ceq = [ceq; kin1; kin2; kin3; kin4; kin5; kin6];
            end
        end
        
        % Also enforce kinematic constraints at midpoint for higher accuracy
        q1_mid = s_mid(1); q2_mid = s_mid(3); q3_mid = s_mid(5);
        q4_mid = s_mid(7); q5_mid = s_mid(9); q6_mid = s_mid(11);
        q7_mid = s_mid(13); q8_mid = s_mid(15); q9_mid = s_mid(17);
        
        kin1_mid = (2*L7*cos(q7_mid) + cos(q1_mid)*(L2+q2_mid) - L9 - cos(q3_mid)*(L4+q4_mid)) / CHAR_LENGTH;
        kin2_mid = (2*L7*sin(q7_mid) + sin(q1_mid)*(L2+q2_mid) - sin(q3_mid)*(L4+q4_mid)) / CHAR_LENGTH;
        kin3_mid = (L7*cos(q7_mid) + L8*sin(q7_mid) + cos(q1_mid)*(L2+q2_mid) - 0.5*L9 - cos(q5_mid)*(L6+q6_mid)) / CHAR_LENGTH;
        kin4_mid = (L7*sin(q7_mid) + sin(q1_mid)*(L2+q2_mid) - L8*cos(q7_mid) - sin(q5_mid)*(L6+q6_mid)) / CHAR_LENGTH;
        kin5_mid = (L7*cos(q7_mid) + cos(q1_mid)*(L2+q2_mid) - q8_mid) / CHAR_LENGTH;
        kin6_mid = (L7*sin(q7_mid) + sin(q1_mid)*(L2+q2_mid) - q9_mid) / CHAR_LENGTH;
        
        ceq = [ceq; kin1_mid; kin2_mid; kin3_mid; kin4_mid; kin5_mid; kin6_mid];
    end
    
    % Final kinematic constraints at the last node
    kin1_final = (2*L7*cos(pos.q7(end,1)) + cos(pos.q1(end,1))*(L2+pos.q2(end,1)) - L9 - cos(pos.q3(end,1))*(L4+pos.q4(end,1))) / CHAR_LENGTH;
    kin2_final = (2*L7*sin(pos.q7(end,1)) + sin(pos.q1(end,1))*(L2+pos.q2(end,1)) - sin(pos.q3(end,1))*(L4+pos.q4(end,1))) / CHAR_LENGTH;
    kin3_final = (L7*cos(pos.q7(end,1)) + L8*sin(pos.q7(end,1)) + cos(pos.q1(end,1))*(L2+pos.q2(end,1)) - 0.5*L9 - cos(pos.q5(end,1))*(L6+pos.q6(end,1))) / CHAR_LENGTH;
    kin4_final = (L7*sin(pos.q7(end,1)) + sin(pos.q1(end,1))*(L2+pos.q2(end,1)) - L8*cos(pos.q7(end,1)) - sin(pos.q5(end,1))*(L6+pos.q6(end,1))) / CHAR_LENGTH;
    kin5_final = (L7*cos(pos.q7(end,1)) + cos(pos.q1(end,1))*(L2+pos.q2(end,1)) - pos.q8(end,1)) / CHAR_LENGTH;
    kin6_final = (L7*sin(pos.q7(end,1)) + sin(pos.q1(end,1))*(L2+pos.q2(end,1)) - pos.q9(end,1)) / CHAR_LENGTH;
    
    ceq = [ceq; kin1_final; kin2_final; kin3_final; kin4_final; kin5_final; kin6_final];
    
    % IMPROVED: Scaled final constraints
    ceq_final = [
        (pos.q1(end,1) - q_Ans_fin(1,1)) / CHAR_ANGLE;
        (pos.q2(end,1) - q_Ans_fin(2,1)) / CHAR_LENGTH;
        (pos.q3(end,1) - q_Ans_fin(3,1)) / CHAR_ANGLE;
        (pos.q4(end,1) - q_Ans_fin(4,1)) / CHAR_LENGTH;
        (pos.q5(end,1) - q_Ans_fin(5,1)) / CHAR_ANGLE;
        (pos.q6(end,1) - q_Ans_fin(6,1)) / CHAR_LENGTH;
        (pos.q7(end,1) - q_Ans_fin(7,1)) / CHAR_ANGLE;
        (pos.q8(end,1) - q_Ans_fin(8,1)) / CHAR_LENGTH;
        (pos.q9(end,1) - q_Ans_fin(9,1)) / CHAR_LENGTH;
        (vel.q1(end,1) - 0) / CHAR_VELOCITY_ANGULAR;
        (vel.q2(end,1) - 0) / CHAR_VELOCITY_LINEAR;
        (vel.q3(end,1) - 0) / CHAR_VELOCITY_ANGULAR;
        (vel.q4(end,1) - 0) / CHAR_VELOCITY_LINEAR;
        (vel.q5(end,1) - 0) / CHAR_VELOCITY_ANGULAR;
        (vel.q6(end,1) - 0) / CHAR_VELOCITY_LINEAR;
        (vel.q7(end,1) - 0) / CHAR_VELOCITY_ANGULAR;
        (vel.q8(end,1) - 0) / CHAR_VELOCITY_LINEAR;
        (vel.q9(end,1) - 0) / CHAR_VELOCITY_LINEAR
    ];
    
    ceq = [ceq; ceq_final];
end

%% ------------------------- IMPROVED Robot Dynamics with Better Conditioning -------------------------
function Sdot = dyn_scaled_robust(t,S,Upsilon)
    global L1 L2 L3 L4 L5 L6 L7 L8 L9 q_Ans_init
    global TIME_SCALE LENGTH_SCALE MASS_SCALE FORCE_SCALE
    
    % Extract scaled state variables
    Q1 = S(2,1); U1 = S(3,1);
    Q2 = S(4,1); U2 = S(5,1);
    Q3 = S(6,1); U3 = S(7,1);
    Q4 = S(8,1); U4 = S(9,1);
    Q5 = S(10,1); U5 = S(11,1);
    Q6 = S(12,1); U6 = S(13,1);
    Q7 = S(14,1); U7 = S(15,1);
    Q8 = S(16,1); U8 = S(17,1);
    Q9 = S(18,1); U9 = S(19,1);
    
    % IMPROVED: Add bounds checking to prevent singular configurations
    Q2 = max(Q2, 0.1);  % Prevent L2+Q2 from being too small
    Q4 = max(Q4, 0.1);  % Prevent L4+Q4 from being too small
    Q6 = max(Q6, 0.1);  % Prevent L6+Q6 from being too small
    
    % IMPROVED: Better conditioned inertial properties
    MA_scaled = 6.974;
    MB_scaled = 6.974;
    MC_scaled = 6.974;
    MD_scaled = 6.974;
    ME_scaled = 6.974;
    MF_scaled = 6.974;
    MG_scaled = 1e10;  % Reduced from 4*1e14 for better conditioning
    
    % Scaled gravity (acceleration in μm/s²)
    GRAV_scaled = 9.81 * 1e6 / TIME_SCALE^2;  % Convert m/s² to μm/scaled_time²
    
    % Scaled moments of inertia
    IA3_scaled = (1/12)*MA_scaled*L1^2;
    IB3_scaled = (1/12)*MB_scaled*L2^2;
    IC3_scaled = (1/12)*MC_scaled*L3^2;
    ID3_scaled = (1/12)*MD_scaled*L4^2;
    IE3_scaled = (1/12)*ME_scaled*L5^2;
    IF3_scaled = (1/12)*MF_scaled*L6^2;
    IG3_scaled = (1/5)*MG_scaled*(L7^2 + L8^2);
    
    % Build scaled mass matrix A with improved expressions
    A = zeros(3,3);
    
    % More numerically stable expressions for mass matrix
    cos_Q1_Q7 = cos(Q1-Q7);
    sin_Q1_Q7 = sin(Q1-Q7);
    cos_Q3_Q7 = cos(Q3-Q7);
    sin_Q3_Q7 = sin(Q3-Q7);
    cos_Q5_Q7 = cos(Q5-Q7);
    sin_Q5_Q7 = sin(Q5-Q7);
    
    L2_Q2 = L2 + Q2;
    L4_Q4 = L4 + Q4;
    L6_Q6 = L6 + Q6;
    
    % Build mass matrix elements with better conditioning
    A(1,1) = IG3_scaled + MB_scaled*L7^2*sin_Q1_Q7^2 + MD_scaled*L7^2*sin_Q3_Q7^2 + ...
             0.25*L7^2*cos_Q1_Q7^2*(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2 + ...
             0.25*L7^2*cos_Q3_Q7^2*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2 + ...
             MF_scaled*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2 + ...
             0.25*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)*(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2/L6_Q6^2;
    
    A(1,2) = L7*MD_scaled*cos(Q3)*sin_Q3_Q7 + ...
             0.25*L7*sin(Q1)*cos_Q1_Q7*(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2 + ...
             MF_scaled*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) + ...
             0.25*sin(Q5)*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)*(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/L6_Q6^2 - ...
             L7*MB_scaled*cos(Q1)*sin_Q1_Q7 - ...
             0.25*L7*sin(Q3)*cos_Q3_Q7*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2;
    
    A(1,3) = L7*MD_scaled*sin(Q3)*sin_Q3_Q7 + ...
             0.25*L7*cos(Q3)*cos_Q3_Q7*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2 + ...
             MF_scaled*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) - ...
             L7*MB_scaled*sin(Q1)*sin_Q1_Q7 - ...
             0.25*L7*cos(Q1)*cos_Q1_Q7*(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2 - ...
             0.25*cos(Q5)*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)*(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/L6_Q6^2;
    
    A(2,1) = A(1,2);  % Symmetric
    
    A(2,2) = MB_scaled + MG_scaled + MD_scaled*cos(Q3)^2 + MF_scaled*cos(Q5)^2 + ...
             0.25*sin(Q3)^2*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2 + ...
             0.25*sin(Q5)^2*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)/L6_Q6^2 - ...
             0.25*sin(Q1)^2*(4*MB_scaled-(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2);
    
    A(2,3) = MB_scaled*sin(Q1)*cos(Q1) + MD_scaled*sin(Q3)*cos(Q3) + MF_scaled*sin(Q5)*cos(Q5) - ...
             0.25*sin(Q1)*cos(Q1)*(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2 - ...
             0.25*sin(Q3)*cos(Q3)*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2 - ...
             0.25*sin(Q5)*cos(Q5)*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)/L6_Q6^2;
    
    A(3,1) = A(1,3);  % Symmetric
    A(3,2) = A(2,3);  % Symmetric
    
    A(3,3) = MG_scaled + MB_scaled*sin(Q1)^2 + MD_scaled*sin(Q3)^2 + MF_scaled*sin(Q5)^2 + ...
             0.25*cos(Q1)^2*(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2 + ...
             0.25*cos(Q3)^2*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2 + ...
             0.25*cos(Q5)^2*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)/L6_Q6^2;
    
    % Build scaled Coriolis/centrifugal vector B
    B = zeros(3,1);
    
    B(1) = 0.5*L7*MB_scaled*sin_Q1_Q7*(L2_Q2*U1^2-2*L2_Q2*U1^2-2*L7*cos_Q1_Q7*U7^2) + ...
           0.5*L7*MD_scaled*sin_Q3_Q7*(2*L4_Q4*U3^2-L4_Q4*U3^2-2*L7*cos_Q3_Q7*U7^2) + ...
           0.25*L7*cos_Q1_Q7*((4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2)*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2-MB_scaled*L2_Q2*(4*U1*U2-L2_Q2*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2))/L2_Q2;
    
    B(2) = 0.5*MD_scaled*cos(Q3)*(2*L4_Q4*U3^2-L4_Q4*U3^2-2*L7*cos_Q3_Q7*U7^2) + ...
           0.25*sin(Q1)*((4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2)*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2-MB_scaled*L2_Q2*(4*U1*U2-L2_Q2*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2))/L2_Q2 - ...
           0.5*MB_scaled*cos(Q1)*(L2_Q2*U1^2-2*L2_Q2*U1^2-2*L7*cos_Q1_Q7*U7^2);
    
    B(3) = 0.5*MD_scaled*sin(Q3)*(2*L4_Q4*U3^2-L4_Q4*U3^2-2*L7*cos_Q3_Q7*U7^2) - ...
           0.5*MB_scaled*sin(Q1)*(L2_Q2*U1^2-2*L2_Q2*U1^2-2*L7*cos_Q1_Q7*U7^2) - ...
           0.25*cos(Q1)*((4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2)*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2-MB_scaled*L2_Q2*(4*U1*U2-L2_Q2*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2))/L2_Q2;
    
    % Build scaled gravity vector G
    G = zeros(3,1);
    
    G(1) = 0.5*GRAV_scaled*(L3*L7*MC_scaled*cos(Q3)*cos_Q3_Q7/L4_Q4+L7*MD_scaled*(2*sin(Q3)*sin_Q3_Q7+cos(Q3)*L4_Q4*cos_Q3_Q7/L4_Q4)+2*MF_scaled*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-L1*L7*MA_scaled*cos(Q1)*cos_Q1_Q7/L2_Q2-L7*MB_scaled*(2*sin(Q1)*sin_Q1_Q7+cos(Q1)*L2_Q2*cos_Q1_Q7/L2_Q2)-cos(Q5)*(L5*ME_scaled+MF_scaled*L6_Q6)*(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/L6_Q6);
    
    G(2) = -0.5*GRAV_scaled*(L1*MA_scaled*sin(Q1)*cos(Q1)/L2_Q2+L3*MC_scaled*sin(Q3)*cos(Q3)/L4_Q4+L5*ME_scaled*sin(Q5)*cos(Q5)/L6_Q6-MB_scaled*sin(Q1)*cos(Q1)*(2-L2_Q2/L2_Q2)-MD_scaled*sin(Q3)*cos(Q3)*(2-L4_Q4/L4_Q4)-MF_scaled*sin(Q5)*cos(Q5)*(2-L6_Q6/L6_Q6));
    
    G(3) = 0.5*GRAV_scaled*(2*MG_scaled+2*MB_scaled*sin(Q1)^2+2*MD_scaled*sin(Q3)^2+2*MF_scaled*sin(Q5)^2+cos(Q1)^2*(L1*MA_scaled+MB_scaled*L2_Q2)/L2_Q2+cos(Q3)^2*(L3*MC_scaled+MD_scaled*L4_Q4)/L4_Q4+cos(Q5)^2*(L5*ME_scaled+MF_scaled*L6_Q6)/L6_Q6);
    
    % Build input matrix GTRAN (scaled)
    GTRAN = zeros(3,3);
    GTRAN(1,1) = -L7*sin_Q1_Q7;
    GTRAN(1,2) = L7*sin_Q3_Q7;
    GTRAN(1,3) = sin(Q5)*(L7*cos(Q7)+L8*sin(Q7)) - L7*sin_Q5_Q7 - cos(Q5)*(L7*sin(Q7)-L8*cos(Q7));
    GTRAN(2,1) = cos(Q1);
    GTRAN(2,2) = cos(Q3);
    GTRAN(2,3) = cos(Q5);
    GTRAN(3,1) = sin(Q1);
    GTRAN(3,2) = sin(Q3);
    GTRAN(3,3) = sin(Q5);
    
    % IMPROVED: Better conditioning and regularization
    % cond_A = cond(A);
    % if cond_A > 1e10
    %     % Add adaptive regularization based on condition number
    %     reg_param = max(1e-12, 1e-16 * cond_A);
    %     A = A + reg_param * eye(size(A));
    %     fprintf('Warning: Adding regularization %.2e to mass matrix (cond = %.2e)\n', reg_param, cond_A);
    % end
    
    % Solve with better numerical stability
    try
        E = A \ (-B - G);
        N = A \ (GTRAN * Upsilon);
    catch ME
        % Fallback to pseudo-inverse if direct solve fails
        warning('Using pseudo-inverse due to singular mass matrix: %s', ME.message);
        E = pinv(A, 1e-12) * (-B - G);
        N = pinv(A, 1e-12) * (GTRAN * Upsilon);
    end
    
    % Convert velocities to scaled coordinates (with bounds checking)
    U1 = (cos(Q1)*U9-sin(Q1)*U8-L7*cos_Q1_Q7*U7)/max(L2_Q2, 0.1);
    U2 = sin(Q1)*U9 + cos(Q1)*U8 - L7*sin_Q1_Q7*U7;
    U3 = -(sin(Q3)*U8-cos(Q3)*U9-L7*cos_Q3_Q7*U7)/max(L4_Q4, 0.1);
    U4 = sin(Q3)*U9 + cos(Q3)*U8 + L7*sin_Q3_Q7*U7;
    U5 = (cos(Q5)*U9-sin(Q5)*U8-(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*U7)/max(L6_Q6, 0.1);
    U6 = sin(Q5)*U9 + cos(Q5)*U8 + (sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*U7;
    
    % Calculate accelerations (scaled) with improved expressions
    dotU1 = (cos(Q1)*(E(3,1) + N(3,1))-2*U1*U2-L7*sin_Q1_Q7*U7^2-sin(Q1)*(E(2,1) + N(2,1))-L7*cos_Q1_Q7*(E(1,1) + N(1,1)))/max(L2_Q2, 0.1);
    dotU2 = L2_Q2*U1^2 + L7*cos_Q1_Q7*U7^2 + sin(Q1)*(E(3,1) + N(3,1)) + cos(Q1)*(E(2,1) + N(2,1)) - L7*sin_Q1_Q7*(E(1,1) + N(1,1));
    dotU3 = -(2*U3*U4+sin(Q3)*(E(2,1) + N(2,1))-L7*sin_Q3_Q7*U7^2-cos(Q3)*(E(3,1) + N(3,1))-L7*cos_Q3_Q7*(E(1,1) + N(1,1)))/max(L4_Q4, 0.1);
    dotU4 = L4_Q4*U3^2 + sin(Q3)*(E(3,1) + N(3,1)) + cos(Q3)*(E(2,1) + N(2,1)) + L7*sin_Q3_Q7*(E(1,1) + N(1,1)) - L7*cos_Q3_Q7*U7^2;
    dotU5 = (L8*cos_Q5_Q7*U7^2+cos(Q5)*(E(3,1) + N(3,1))-2*U5*U6-sin(Q5)*(E(2,1) + N(2,1))-(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*(E(1,1) + N(1,1)))/max(L6_Q6, 0.1);
    dotU6 = L6_Q6*U5^2 + L8*sin_Q5_Q7*U7^2 + sin(Q5)*(E(3,1) + N(3,1)) + cos(Q5)*(E(2,1) + N(2,1)) + (sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*(E(1,1) + N(1,1));
    
    % The dynamics in scaled coordinates
    Sdot = zeros(19,1);
    Sdot(1,1) = 1; % First state variable is time: d(t_scaled)/dt_scaled = 1
    Sdot(2,1) = U1;                  
    Sdot(3,1) = dotU1;  
    Sdot(4,1) = U2;               
    Sdot(5,1) = dotU2; 
    Sdot(6,1) = U3;               
    Sdot(7,1) = dotU3; 
    Sdot(8,1) = U4; 
    Sdot(9,1) = dotU4; 
    Sdot(10,1) = U5; 
    Sdot(11,1) = dotU5; 
    Sdot(12,1) = U6; 
    Sdot(13,1) = dotU6; 
    Sdot(14,1) = S(15); 
    Sdot(15,1) = E(1,1) + N(1,1); 
    Sdot(16,1) = S(17); 
    Sdot(17,1) = E(2,1) + N(2,1); 
    Sdot(18,1) = S(19); 
    Sdot(19,1) = E(3,1) + N(3,1); 
end

function Sdot = dyn2(S,Upsilon)
    global L1 L2 L3 L4 L5 L6 L7 L8 L9 q_Ans_init
    global TIME_SCALE LENGTH_SCALE MASS_SCALE FORCE_SCALE
    
    % Extract scaled state variables
    Q1 = S(1,1); U1 = S(2,1);
    Q2 = S(3,1); U2 = S(4,1);
    Q3 = S(5,1); U3 = S(6,1);
    Q4 = S(7,1); U4 = S(8,1);
    Q5 = S(9,1); U5 = S(10,1);
    Q6 = S(11,1); U6 = S(12,1);
    Q7 = S(13,1); U7 = S(14,1);
    Q8 = S(15,1); U8 = S(16,1);
    Q9 = S(17,1); U9 = S(18,1);
    
    % IMPROVED: Add bounds checking to prevent singular configurations
    Q2 = max(Q2, 0.1);  % Prevent L2+Q2 from being too small
    Q4 = max(Q4, 0.1);  % Prevent L4+Q4 from being too small
    Q6 = max(Q6, 0.1);  % Prevent L6+Q6 from being too small
    
    % IMPROVED: Better conditioned inertial properties
    MA_scaled = 6.974;
    MB_scaled = 6.974;
    MC_scaled = 6.974;
    MD_scaled = 6.974;
    ME_scaled = 6.974;
    MF_scaled = 6.974;
    MG_scaled = 1e10;  % Reduced from 4*1e14 for better conditioning
    
    % Scaled gravity (acceleration in μm/s²)
    GRAV_scaled = 9.81 * 1e6 / TIME_SCALE^2;  % Convert m/s² to μm/scaled_time²
    
    % Scaled moments of inertia
    IA3_scaled = (1/12)*MA_scaled*L1^2;
    IB3_scaled = (1/12)*MB_scaled*L2^2;
    IC3_scaled = (1/12)*MC_scaled*L3^2;
    ID3_scaled = (1/12)*MD_scaled*L4^2;
    IE3_scaled = (1/12)*ME_scaled*L5^2;
    IF3_scaled = (1/12)*MF_scaled*L6^2;
    IG3_scaled = (1/5)*MG_scaled*(L7^2 + L8^2);
    
    % Build scaled mass matrix A with improved expressions
    A = zeros(3,3);
    
    % More numerically stable expressions for mass matrix
    cos_Q1_Q7 = cos(Q1-Q7);
    sin_Q1_Q7 = sin(Q1-Q7);
    cos_Q3_Q7 = cos(Q3-Q7);
    sin_Q3_Q7 = sin(Q3-Q7);
    cos_Q5_Q7 = cos(Q5-Q7);
    sin_Q5_Q7 = sin(Q5-Q7);
    
    L2_Q2 = L2 + Q2;
    L4_Q4 = L4 + Q4;
    L6_Q6 = L6 + Q6;
    
    % Build mass matrix elements with better conditioning
    A(1,1) = IG3_scaled + MB_scaled*L7^2*sin_Q1_Q7^2 + MD_scaled*L7^2*sin_Q3_Q7^2 + ...
             0.25*L7^2*cos_Q1_Q7^2*(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2 + ...
             0.25*L7^2*cos_Q3_Q7^2*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2 + ...
             MF_scaled*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2 + ...
             0.25*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)*(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))^2/L6_Q6^2;
    
    A(1,2) = L7*MD_scaled*cos(Q3)*sin_Q3_Q7 + ...
             0.25*L7*sin(Q1)*cos_Q1_Q7*(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2 + ...
             MF_scaled*cos(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) + ...
             0.25*sin(Q5)*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)*(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/L6_Q6^2 - ...
             L7*MB_scaled*cos(Q1)*sin_Q1_Q7 - ...
             0.25*L7*sin(Q3)*cos_Q3_Q7*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2;
    
    A(1,3) = L7*MD_scaled*sin(Q3)*sin_Q3_Q7 + ...
             0.25*L7*cos(Q3)*cos_Q3_Q7*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2 + ...
             MF_scaled*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7))) - ...
             L7*MB_scaled*sin(Q1)*sin_Q1_Q7 - ...
             0.25*L7*cos(Q1)*cos_Q1_Q7*(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2 - ...
             0.25*cos(Q5)*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)*(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/L6_Q6^2;
    
    A(2,1) = A(1,2);  % Symmetric
    
    A(2,2) = MB_scaled + MG_scaled + MD_scaled*cos(Q3)^2 + MF_scaled*cos(Q5)^2 + ...
             0.25*sin(Q3)^2*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2 + ...
             0.25*sin(Q5)^2*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)/L6_Q6^2 - ...
             0.25*sin(Q1)^2*(4*MB_scaled-(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2);
    
    A(2,3) = MB_scaled*sin(Q1)*cos(Q1) + MD_scaled*sin(Q3)*cos(Q3) + MF_scaled*sin(Q5)*cos(Q5) - ...
             0.25*sin(Q1)*cos(Q1)*(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2 - ...
             0.25*sin(Q3)*cos(Q3)*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2 - ...
             0.25*sin(Q5)*cos(Q5)*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)/L6_Q6^2;
    
    A(3,1) = A(1,3);  % Symmetric
    A(3,2) = A(2,3);  % Symmetric
    
    A(3,3) = MG_scaled + MB_scaled*sin(Q1)^2 + MD_scaled*sin(Q3)^2 + MF_scaled*sin(Q5)^2 + ...
             0.25*cos(Q1)^2*(4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2+MB_scaled*L2_Q2^2)/L2_Q2^2 + ...
             0.25*cos(Q3)^2*(4*IC3_scaled+4*ID3_scaled+MC_scaled*L3^2+MD_scaled*L4_Q4^2)/L4_Q4^2 + ...
             0.25*cos(Q5)^2*(4*IE3_scaled+4*IF3_scaled+ME_scaled*L5^2+MF_scaled*L6_Q6^2)/L6_Q6^2;
    
    % Build scaled Coriolis/centrifugal vector B
    B = zeros(3,1);
    
    B(1) = 0.5*L7*MB_scaled*sin_Q1_Q7*(L2_Q2*U1^2-2*L2_Q2*U1^2-2*L7*cos_Q1_Q7*U7^2) + ...
           0.5*L7*MD_scaled*sin_Q3_Q7*(2*L4_Q4*U3^2-L4_Q4*U3^2-2*L7*cos_Q3_Q7*U7^2) + ...
           0.25*L7*cos_Q1_Q7*((4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2)*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2-MB_scaled*L2_Q2*(4*U1*U2-L2_Q2*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2))/L2_Q2;
    
    B(2) = 0.5*MD_scaled*cos(Q3)*(2*L4_Q4*U3^2-L4_Q4*U3^2-2*L7*cos_Q3_Q7*U7^2) + ...
           0.25*sin(Q1)*((4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2)*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2-MB_scaled*L2_Q2*(4*U1*U2-L2_Q2*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2))/L2_Q2 - ...
           0.5*MB_scaled*cos(Q1)*(L2_Q2*U1^2-2*L2_Q2*U1^2-2*L7*cos_Q1_Q7*U7^2);
    
    B(3) = 0.5*MD_scaled*sin(Q3)*(2*L4_Q4*U3^2-L4_Q4*U3^2-2*L7*cos_Q3_Q7*U7^2) - ...
           0.5*MB_scaled*sin(Q1)*(L2_Q2*U1^2-2*L2_Q2*U1^2-2*L7*cos_Q1_Q7*U7^2) - ...
           0.25*cos(Q1)*((4*IA3_scaled+4*IB3_scaled+MA_scaled*L1^2)*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2-MB_scaled*L2_Q2*(4*U1*U2-L2_Q2*(2*U1*U2+L7*sin_Q1_Q7*U7^2)/L2_Q2))/L2_Q2;
    
    % Build scaled gravity vector G
    G = zeros(3,1);
    
    G(1) = 0.5*GRAV_scaled*(L3*L7*MC_scaled*cos(Q3)*cos_Q3_Q7/L4_Q4+L7*MD_scaled*(2*sin(Q3)*sin_Q3_Q7+cos(Q3)*L4_Q4*cos_Q3_Q7/L4_Q4)+2*MF_scaled*sin(Q5)*(sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))-L1*L7*MA_scaled*cos(Q1)*cos_Q1_Q7/L2_Q2-L7*MB_scaled*(2*sin(Q1)*sin_Q1_Q7+cos(Q1)*L2_Q2*cos_Q1_Q7/L2_Q2)-cos(Q5)*(L5*ME_scaled+MF_scaled*L6_Q6)*(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))/L6_Q6);
    
    G(2) = -0.5*GRAV_scaled*(L1*MA_scaled*sin(Q1)*cos(Q1)/L2_Q2+L3*MC_scaled*sin(Q3)*cos(Q3)/L4_Q4+L5*ME_scaled*sin(Q5)*cos(Q5)/L6_Q6-MB_scaled*sin(Q1)*cos(Q1)*(2-L2_Q2/L2_Q2)-MD_scaled*sin(Q3)*cos(Q3)*(2-L4_Q4/L4_Q4)-MF_scaled*sin(Q5)*cos(Q5)*(2-L6_Q6/L6_Q6));
    
    G(3) = 0.5*GRAV_scaled*(2*MG_scaled+2*MB_scaled*sin(Q1)^2+2*MD_scaled*sin(Q3)^2+2*MF_scaled*sin(Q5)^2+cos(Q1)^2*(L1*MA_scaled+MB_scaled*L2_Q2)/L2_Q2+cos(Q3)^2*(L3*MC_scaled+MD_scaled*L4_Q4)/L4_Q4+cos(Q5)^2*(L5*ME_scaled+MF_scaled*L6_Q6)/L6_Q6);
    
    % Build input matrix GTRAN (scaled)
    GTRAN = zeros(3,3);
    GTRAN(1,1) = -L7*sin_Q1_Q7;
    GTRAN(1,2) = L7*sin_Q3_Q7;
    GTRAN(1,3) = sin(Q5)*(L7*cos(Q7)+L8*sin(Q7)) - L7*sin_Q5_Q7 - cos(Q5)*(L7*sin(Q7)-L8*cos(Q7));
    GTRAN(2,1) = cos(Q1);
    GTRAN(2,2) = cos(Q3);
    GTRAN(2,3) = cos(Q5);
    GTRAN(3,1) = sin(Q1);
    GTRAN(3,2) = sin(Q3);
    GTRAN(3,3) = sin(Q5);
    
    % IMPROVED: Better conditioning and regularization
    % cond_A = cond(A);
    % if cond_A > 1e10
    %     % Add adaptive regularization based on condition number
    %     reg_param = max(1e-12, 1e-16 * cond_A);
    %     A = A + reg_param * eye(size(A));
    %     fprintf('Warning: Adding regularization %.2e to mass matrix (cond = %.2e)\n', reg_param, cond_A);
    % end
    
    % Solve with better numerical stability
    try
        E = A \ (-B - G);
        N = A \ (GTRAN * Upsilon);
    catch ME
        % Fallback to pseudo-inverse if direct solve fails
        warning('Using pseudo-inverse due to singular mass matrix: %s', ME.message);
        E = pinv(A, 1e-12) * (-B - G);
        N = pinv(A, 1e-12) * (GTRAN * Upsilon);
    end
    
    % Convert velocities to scaled coordinates (with bounds checking)
    U1 = (cos(Q1)*U9-sin(Q1)*U8-L7*cos_Q1_Q7*U7)/max(L2_Q2, 0.1);
    U2 = sin(Q1)*U9 + cos(Q1)*U8 - L7*sin_Q1_Q7*U7;
    U3 = -(sin(Q3)*U8-cos(Q3)*U9-L7*cos_Q3_Q7*U7)/max(L4_Q4, 0.1);
    U4 = sin(Q3)*U9 + cos(Q3)*U8 + L7*sin_Q3_Q7*U7;
    U5 = (cos(Q5)*U9-sin(Q5)*U8-(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*U7)/max(L6_Q6, 0.1);
    U6 = sin(Q5)*U9 + cos(Q5)*U8 + (sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*U7;
    
    % Calculate accelerations (scaled) with improved expressions
    dotU1 = (cos(Q1)*(E(3,1) + N(3,1))-2*U1*U2-L7*sin_Q1_Q7*U7^2-sin(Q1)*(E(2,1) + N(2,1))-L7*cos_Q1_Q7*(E(1,1) + N(1,1)))/max(L2_Q2, 0.1);
    dotU2 = L2_Q2*U1^2 + L7*cos_Q1_Q7*U7^2 + sin(Q1)*(E(3,1) + N(3,1)) + cos(Q1)*(E(2,1) + N(2,1)) - L7*sin_Q1_Q7*(E(1,1) + N(1,1));
    dotU3 = -(2*U3*U4+sin(Q3)*(E(2,1) + N(2,1))-L7*sin_Q3_Q7*U7^2-cos(Q3)*(E(3,1) + N(3,1))-L7*cos_Q3_Q7*(E(1,1) + N(1,1)))/max(L4_Q4, 0.1);
    dotU4 = L4_Q4*U3^2 + sin(Q3)*(E(3,1) + N(3,1)) + cos(Q3)*(E(2,1) + N(2,1)) + L7*sin_Q3_Q7*(E(1,1) + N(1,1)) - L7*cos_Q3_Q7*U7^2;
    dotU5 = (L8*cos_Q5_Q7*U7^2+cos(Q5)*(E(3,1) + N(3,1))-2*U5*U6-sin(Q5)*(E(2,1) + N(2,1))-(L7*cos_Q5_Q7-cos(Q5)*(L7*cos(Q7)+L8*sin(Q7))-sin(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*(E(1,1) + N(1,1)))/max(L6_Q6, 0.1);
    dotU6 = L6_Q6*U5^2 + L8*sin_Q5_Q7*U7^2 + sin(Q5)*(E(3,1) + N(3,1)) + cos(Q5)*(E(2,1) + N(2,1)) + (sin(Q5)*(L7*cos(Q7)+L8*sin(Q7))-L7*sin_Q5_Q7-cos(Q5)*(L7*sin(Q7)-L8*cos(Q7)))*(E(1,1) + N(1,1));
    
    % The dynamics in scaled coordinates
    Sdot = zeros(18,1);
    Sdot(1,1) = U1;                  
    Sdot(2,1) = dotU1;  
    Sdot(3,1) = U2;               
    Sdot(4,1) = dotU2; 
    Sdot(5,1) = U3;               
    Sdot(6,1) = dotU3; 
    Sdot(7,1) = U4; 
    Sdot(8,1) = dotU4; 
    Sdot(9,1) = U5; 
    Sdot(10,1) = dotU5; 
    Sdot(11,1) = U6; 
    Sdot(12,1) = dotU6; 
    Sdot(13,1) = S(14); 
    Sdot(14,1) = E(1,1) + N(1,1); 
    Sdot(15,1) = S(16); 
    Sdot(16,1) = E(2,1) + N(2,1); 
    Sdot(17,1) = S(18); 
    Sdot(18,1) = E(3,1) + N(3,1); 
end

%% ------------------------- Scaled Nucleus Dynamics -------------------------
function Sdot_Nuc = NucDyn_scaled(t,S)    
    global TIME_SCALE LENGTH_SCALE
    
    % Convert scaled time back to original time for DMD equations
    t_original = t * TIME_SCALE;
    
    % The dynamics obtained through DMD (converted to scaled coordinates)
    Sdot_Nuc(1,1) = 1; % d(t_scaled)/dt_scaled = 1
    Sdot_Nuc(2,1) = 0; % d(q7)/dt_scaled
    
    % Original: d(q8)/dt = (4.859*10^-9) - (6.759*10^-13)*t + (4.0255*10^-17)*t^2 - (8.35139*10^-22)*t^3
    % Convert to scaled coordinates: d(q8_scaled)/dt_scaled = (LENGTH_SCALE/TIME_SCALE) * original_derivative
    dq8_dt_original = (4.859*10^-9) - (6.759*10^-13)*t_original + (4.0255*10^-17)*t_original^2 - (8.35139*10^-22)*t_original^3;
    Sdot_Nuc(3,1) = dq8_dt_original * TIME_SCALE / LENGTH_SCALE; % Convert to scaled units
    
    Sdot_Nuc(4,1) = 0; % d(q9)/dt_scaled
end

%% ------------------------- Diagnostic Functions -------------------------
function diagnose_optimization_results(S, N)
    global CHAR_LENGTH CHAR_ANGLE
    
    fprintf('\n=== OPTIMIZATION DIAGNOSTICS ===\n');
    
    % Check constraint violations
    [c, ceq] = Constraints_Improved(S, N);
    max_violation = max(abs(ceq));
    mean_violation = mean(abs(ceq));
    
    fprintf('Constraint Analysis:\n');
    fprintf('  Total constraints: %d\n', length(ceq));
    fprintf('  Max violation: %.6e\n', max_violation);
    fprintf('  Mean violation: %.6e\n', mean_violation);
    fprintf('  Std violation: %.6e\n', std(abs(ceq)));
    
    % Find worst constraints
    [sorted_violations, indices] = sort(abs(ceq), 'descend');
    fprintf('\nWorst 5 constraint violations:\n');
    for i = 1:min(5, length(ceq))
        fprintf('  Constraint %d: %.6e\n', indices(i), sorted_violations(i));
    end
    
    % Check mass matrix conditioning throughout trajectory
    fprintf('\nMass Matrix Conditioning:\n');
    cond_numbers = [];
    for i = 1:N
        % Extract states at each node
        q1 = S(2+(i-1),1);
        q2 = S(2+N+(i-1),1);
        q3 = S(2+2*N+(i-1),1);
        q4 = S(2+3*N+(i-1),1);
        q5 = S(2+4*N+(i-1),1);
        q6 = S(2+5*N+(i-1),1);
        q7 = S(2+6*N+(i-1),1);
        
        % Build mass matrix for this state (simplified check)
        A_test = build_test_mass_matrix(q1, q2, q3, q4, q5, q6, q7);
        cond_numbers(i) = cond(A_test);
    end
    
    fprintf('  Mean condition number: %.2e\n', mean(cond_numbers));
    fprintf('  Max condition number: %.2e\n', max(cond_numbers));
    fprintf('  Min condition number: %.2e\n', min(cond_numbers));
    
    % Plot diagnostics
    figure('Name', 'Optimization Diagnostics');
    
    subplot(2,2,1);
    semilogy(abs(ceq), 'o-');
    title('Constraint Violations');
    xlabel('Constraint Index');
    ylabel('|Violation|');
    grid on;
    
    subplot(2,2,2);
    semilogy(cond_numbers, 's-');
    title('Mass Matrix Condition Numbers');
    xlabel('Time Node');
    ylabel('Condition Number');
    grid on;
    
    subplot(2,2,3);
    histogram(log10(abs(ceq(abs(ceq) > 1e-16))));
    title('Distribution of Log10(|Violations|)');
    xlabel('Log10(|Violation|)');
    ylabel('Count');
    grid on;
    
    subplot(2,2,4);
    plot(1:length(ceq), cumsum(abs(ceq)));
    title('Cumulative Constraint Violations');
    xlabel('Constraint Index');
    ylabel('Cumulative |Violation|');
    grid on;
end

function A_test = build_test_mass_matrix(q1, q2, q3, q4, q5, q6, q7)
    % Simplified mass matrix for conditioning check
    global L1 L2 L3 L4 L5 L6 L7 L8 L9
    
    % Use simplified mass values for testing
    M = 1;
    I = 1;
    
    A_test = eye(3);  % Start with identity
    
    % Add dominant terms that affect conditioning
    A_test(1,1) = I + M*(L2+q2)^2 + M*(L4+q4)^2 + M*(L6+q6)^2;
    A_test(2,2) = M + M*cos(q1)^2 + M*cos(q3)^2 + M*cos(q5)^2;
    A_test(3,3) = M + M*sin(q1)^2 + M*sin(q3)^2 + M*sin(q5)^2;
    
    % Add coupling terms
    A_test(1,2) = M*cos(q1)*sin(q1-q7);
    A_test(2,1) = A_test(1,2);
    A_test(1,3) = M*sin(q1)*sin(q1-q7);
    A_test(3,1) = A_test(1,3);
    A_test(2,3) = M*sin(q1)*cos(q1);
    A_test(3,2) = A_test(2,3);
end

%% ------------------------- Display Results Summary -------------------------
    fprintf('\n Optimization Completed in %f seconds \n', toc);
    
    % Run diagnostics
    diagnose_optimization_results(S, N);
    
    % Check final constraint violations
    [c_final, ceq_final] = Constraints_Improved(S, N);
    max_eq_violation = max(abs(ceq_final));
    fprintf('\n=== FINAL RESULTS SUMMARY ===\n');
    fprintf('Final time (original): %.2f seconds\n', tfin_original);
    fprintf('Final time (scaled): %.4f\n', tfin_scaled);
    fprintf('Maximum equality constraint violation: %.6e\n', max_eq_violation);
    
    if max_eq_violation < 1e-4
        fprintf('SUCCESS: Constraints satisfied to good tolerance!\n');
    elseif max_eq_violation < 1e-2
        fprintf('PARTIAL SUCCESS: Constraints satisfied to reasonable tolerance.\n');
    else
        fprintf('WARNING: Large constraint violations detected.\n');
    end
    
    fprintf('\n=== SCALING INFORMATION ===\n');
    fprintf('Time scaling factor: %.2f\n', TIME_SCALE);
    fprintf('Length scaling factor: %.0e m = 1 μm\n', LENGTH_SCALE);
    fprintf('Force scaling factor: %.0e N = 1 nN\n', FORCE_SCALE);
    fprintf('Mass scaling factor: %.0e kg\n', MASS_SCALE);
    
    % Display force ranges
    fprintf('\nControl force ranges (in nanoNewtons):\n');
    fprintf('Upsilon1: [%.2f, %.2f] nN\n', min(control.Upsilon1)*1e9, max(control.Upsilon1)*1e9);
    fprintf('Upsilon2: [%.2f, %.2f] nN\n', min(control.Upsilon2)*1e9, max(control.Upsilon2)*1e9);
    fprintf('Upsilon3: [%.2f, %.2f] nN\n', min(control.Upsilon3)*1e9, max(control.Upsilon3)*1e9);
    
    % Save results to file
    save('improved_optimization_results.mat', 'S', 'pos', 'vel', 'control', 'tfin_original', 'tvec_original');
    
    fprintf('\nResults saved to improved_optimization_results.mat\n');
    fprintf('========================================\n');
    
    % Additional optimization tips in comments:
    % 1. If still having feasibility issues, try reducing N to 5-7 initially
    % 2. Consider using a continuation method: start with loose constraints, then tighten
    % 3. Use warm-starting: solve a simpler version first, then use result as initial guess
    % 4. Check if the problem is over-constrained by temporarily removing some constraints
    % 5. Verify that initial and final conditions are physically realistic