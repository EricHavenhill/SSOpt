%% |==============================================================================================|%
  %|      Filename: StressFiberOptimizationManyAdhesionsV10.m                                     |%
  %|                                                                                              |%
  %|      Author  : Eric Havenhill, Ph.D. (Modified - Fixed Negative Areas)                       |%
  %|                Faculty Affiliate                                                             |%
  %|                Colorado State University                                                     |%
  %|                                                                                              |%
  %|      Purpose : This code performs a structural optimization (sizing) on a cell with 7        |%
  %                 elements serving as stress fibers                                             |%
  %|                                                                                              |%
  %|      Notes   : * FIXED: Added proper bounds to prevent negative areas                        |%
  %|                * FIXED: Improved constraint formulation                                      |%
  %|                * FIXED: Better initial conditions and scaling                                |%
  %|                * FIXED: Anonymous function call to MyNLObj                                   |%
  %|                * MODIFIED: Bell curve initial guess for areas                                |%
  %|                * MODIFIED: Peak set to 4.5 μm²                                               |%
  %|                * V9: Fixed finite difference warnings with better bounds and scaling         |%
  %|                * V9: Fixed infeasibility issues with realistic parameters                    |%
%% |==============================================================================================|%

%% =========================================== Setup ==============================================

    clc; clear all; close all;

    format long
        
    global E Fa sigma_max L_mid theta L Ne Nn
    
    Ne = 7;
    Nn = Ne +1;
    
    % Small Elasticity
    %E = 0.5*10^6;  % Use softer fibers for protruding case
    
    % Mid. Elasticity 
    %E = 1.45*10^6;
    
    % Large Elasticity
    E = 2*10^6;
    
    % FIXED: More realistic cellular force magnitudes
    %Fa = 10*10^-9;  % 10 nN - more typical for cellular mechanics
    Fa = 1*10^-9;
    %Fa = 0;

    sigma_max = 10*10^6;
    
    % For the non-protruding (symmetric) case:
    %theta = (pi/180)*linspace(15,180-15,Ne)';

    % For the protruding (offset) case:
    theta = (pi/180)*linspace(90,180-15,Ne)';
    
    L_mid = 26*10^-6;
    
    for n = 1 : Ne
        L(1,n) = L_mid/sin(theta(n,1));
    end
    
    % MODIFIED: Bell curve initial guess with peak at 4.5 μm²
    % Create bell curve centered at the 4th element (middle fiber)
    center_index = 4;  % Central fiber (4th out of 7)
    peak_area = 4.5e-12;  % Peak area at center (4.5 μm² = 4.5e-12 m²)
    min_area = 0.5e-12;  % Minimum area at ends
    
    % Create indices for bell curve calculation
    indices = 1:Ne;
    
    % Bell curve parameters
    sigma_bell = 2;  % Controls width of bell curve
    
    % Calculate bell curve values
    bell_values = exp(-0.5 * ((indices - center_index) / sigma_bell).^2);
    
    % Scale bell curve to desired range
    A_init = min_area + (peak_area - min_area) * bell_values';
    
    % Display the bell curve initial areas
    fprintf('Bell curve initial areas (m²):\n');
    for i = 1:Ne
        fprintf('Fiber %d: %e m² (%.4f μm²)\n', i, A_init(i), A_init(i)*1e12);
    end
    fprintf('\n');
    
    u_init_X = 2*10^-6;    
    u_init_Y = -1*10^-6;
    u_init = [u_init_X; u_init_Y];  

    X0 = [A_init; u_init.*ones(2,1)]; 

%% ======================================== Optimization ===========================================

    % FIXED: More relaxed bounds for protruding case
    % Lower bounds: areas must be positive with reasonable minimum
    lb = [1e-14*ones(Ne,1);      % Minimum area bound (very small but positive)
          0e-6;                 % ux lower bound (allow negative)
          -20e-6];                % uy lower bound
    
    % Upper bounds: reasonable maximums with more room
    ub = [20e-12*ones(Ne,1);     % Maximum area bound (20 μm²)
          20e-6;                  % ux upper bound (allow positive)
          0e-6];                 % uy upper bound (allow positive - more flexible)
    
    % --------------------------------- Optimization with fmincon ---------------------------------
    
    % Suppress finite difference warnings (they're informational, not errors)
    warning('off', 'optim:fmincon:SwitchingToMediumScale');
    
    options = optimoptions(@fmincon,... 
        'TolFun',1E-6, ...              
        'TolCon',1E-6,...               
        'TolX',1E-8,...                 
        'MaxFunEvals',50000, ...        
        'MaxIter', 2000, ...
        'Display', 'iter', ...
        'Algorithm', 'interior-point',...  % Changed from 'sqp' - handles bounds better
        'ScaleProblem', true,...
        'FiniteDifferenceType', 'central',...  % More accurate derivatives
        'OptimalityTolerance', 1e-6,...
        'ConstraintTolerance', 1e-6,...
        'StepTolerance', 1e-10);
        % Removed FiniteDifferenceStepSize to use automatic scaling
    
    % % FIXED: Include bounds in fmincon call
    % [S1,fval,exitflag,output] = fmincon(@MyNLObj, X0, [], [], [], [], lb, ub, @MyNLCon, options);
    % 
    % % Turn warnings back on
    warning('on', 'all');

    % ------------------------------- Optimization with MultiStart -------------------------------
    
    % FIXED: Corrected anonymous function call - pass the argument to MyNLObj
    problem = createOptimProblem('fmincon', 'x0', X0, 'objective', @(S) MyNLObj(S), ...
                                    'nonlcon', @MyNLCon, 'lb', lb, 'ub', ub, 'options', options);

    ms = MultiStart;
    [S1, fval, exitflag, output, solutions] = run(ms, problem, 10); % 100 random starts

    % ------------------------------- Optimization with GlobalSearch -------------------------------

    % gs = GlobalSearch('NumTrialPoints', 1000, 'NumStageOnePoints', 200);
    % [S1, fval, exitflag, output] = run(gs, problem);
    
    % ------------------------- Display Results -------------------------
    % Display convergence information
    fprintf('\n=============== CONVERGENCE INFO ===============\n');
    fprintf('Exit flag: %d\n', exitflag);
    fprintf('Function evaluations: %d\n', output.funcCount);
    % fprintf('Iterations: %d\n', output.iterations);
    
    % Check exit flag meaning
    if exitflag > 0
        fprintf('Status: Successfully converged\n');
    elseif exitflag == 0
        fprintf('Status: Maximum iterations/evaluations reached\n');
    else
        fprintf('Status: Did not converge (possibly infeasible)\n');
    end
    
    % Verify equilibrium is satisfied
    [g, geq] = MyNLCon(S1);
    fprintf('Max inequality constraint violation: %e\n', max(max(g), 0));
    fprintf('Equality constraint violations: [%e, %e]\n', geq(1), geq(2));
    
    for i = 1: Ne
        strain(i) = (cos(theta(i))*S1(end-1) + sin(theta(i))*S1(end))/L(i);
        FiberForce(i) = S1(i)*E*strain(i);
        FiberStress(i) = E*strain(i);
    end

    fprintf('\n=============== OPTIMIZATION RESULTS ===============\n');
    fprintf('Optimization completed with objective value: %e\n', fval);
    fprintf('\n');
    fprintf('Optimized areas (μm²):\n');
    for i = 1:Ne
        fprintf('  Fiber %d: %.6f μm²\n', i, S1(i)*1e12);
    end
    fprintf('\n');
    fprintf('Total volume: %e m³\n', sum(S1(1:Ne).*L'));
    fprintf('\n');
    fprintf('Optimized displacements:\n');
    fprintf('  ux = %e m (%.4f μm)\n', S1(end-1), S1(end-1)*1e6);
    fprintf('  uy = %e m (%.4f μm)\n', S1(end), S1(end)*1e6);
    fprintf('\n');
    fprintf('Fiber strains:\n');
    for i = 1:Ne
        fprintf('  Fiber %d: %e (%.6f%%)\n', i, strain(i), strain(i)*100);
    end
    fprintf('\n');
    fprintf('Fiber forces (N):\n');
    for i = 1:Ne
        fprintf('  Fiber %d: %e N\n', i, FiberForce(i));
    end
    fprintf('Total fiber force (horizontal): %e N\n', sum(FiberForce.*cos(theta')));
    fprintf('Total fiber force (vertical): %e N\n', sum(FiberForce.*sin(theta')));
    fprintf('\n');
    fprintf('Fiber stresses (Pa):\n');
    for i = 1:Ne
        fprintf('  Fiber %d: %e Pa\n', i, FiberStress(i));
    end
    
    % FIXED: Check for negative areas and warn
    if any(S1(1:Ne) < 0)
        warning('Negative areas detected! This should not happen with proper bounds.');
        fprintf('Negative areas at indices: %s\n', num2str(find(S1(1:Ne) < 0)'));
    end
    
    % Check if any areas are at bounds
    tol = 1e-13;
    at_lower = S1(1:Ne) - lb(1:Ne) < tol;
    at_upper = ub(1:Ne) - S1(1:Ne) < tol;
    if any(at_lower)
        fprintf('\nWarning: Fibers at lower bound: %s\n', num2str(find(at_lower)'));
    end
    if any(at_upper)
        fprintf('\nWarning: Fibers at upper bound: %s\n', num2str(find(at_upper)'));
    end

%% ======================================= Cost Function ==========================================
    
    function J = MyNLObj(X)
        global E Fa sigma_max L_mid theta L Ne Nn
        
        % Minimize total volume (area * length)
        J = 1E0*L*X(1:Ne,1);
    end

%% ==================================== Constraints ====================================
    
    function [g,geq] = MyNLCon(X)
        global E Fa sigma_max L_mid theta L Ne Nn
       
        % Initialize constraint arrays
        g1 = zeros(Ne,1);
        g2 = zeros(Ne,1);
        g3 = zeros(Ne,1);
        
        for i=1:Ne
            % Calculate strain for this fiber
            strain = (cos(theta(i))*X(end-1,1) + sin(theta(i))*X(end,1))/L(i);
            stress = E * strain;
            
            % Stress constraints (tension and compression limits)
            g1(i,1) = stress - sigma_max;      % Stress <= sigma_max
            g2(i,1) = -stress - sigma_max;     % Stress >= -sigma_max
            
            % FIXED: Buckling constraint - only active under compression
            if strain < 0  % Compression
                % Critical buckling stress for circular cross-section
                sigma_critical = pi * E * X(i,1) / (4 * L(i)^2);
                g3(i,1) = -stress - sigma_critical;  % |compression stress| <= critical stress
            else
                g3(i,1) = -1;  % Always satisfied in tension
            end
        end
        
        % FIXED: Area constraints now handled by bounds, but keep as backup
        % (these should be inactive due to bounds)
        g4 = -X(1:Ne,1) + 1e-15;  % Areas >= 1e-15 (essentially zero)
        
        % Combine all inequality constraints
        g = [g1; g2; g3; g4];
        
        % FIXED: Equilibrium constraints from FEM formulation
        k11 = 0;
        k12 = 0;
        k22 = 0;
        
        for i = 1:Ne
            k11 = k11 + X(i,1)*E*cos(theta(i))^2/L(i);
            k12 = k12 + X(i,1)*E*cos(theta(i))*sin(theta(i))/L(i);
            k22 = k22 + X(i,1)*E*sin(theta(i))^2/L(i);
        end
        
        k21 = k12;  % Symmetry
        
        % FIXED: Force balance equations
        g_grav = 9.81;
        m_cell = 2.290*10^-11;

        geq1 = Fa - k11*X(end-1,1) - k12*X(end,1);     % Horizontal equilibrium
        geq2 = -m_cell*g_grav - k21*X(end-1,1) - k22*X(end,1);  % Vertical equilibrium
        
        geq = [geq1; geq2];
        
    end