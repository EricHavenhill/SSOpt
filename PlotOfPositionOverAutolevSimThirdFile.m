% Run this one third?
someName = load( 'cellsimplifiedmodv4_autolev_2_matlab.1' );
% someName = load( 'cellsimplifiedmodv4druggedm_autolev_2_matlab.1' );

figure(1)
%grid on;
    % Set figure dimensions (width x height in pixels)
    fig_width = 600;
    fig_height = 500;
    set(gcf, 'Position', [100, 100, fig_width, fig_height]);
    % Set font size (customizable)
    font_size = 16.5; % Change this value as needed
plot( someName(:,1), someName(:,8), 'b', someName(:,1), someName(:,9), 'r',someName(:,1), someName(:,10), 'k' )
xlabel('t (s)'); ylabel('Actuator Forces', 'Interpreter','latex');
%grid on;
% legend('$\Upsilon_{1}$','$\Upsilon_{2}$','$\Upsilon_{3}$','Interpreter','latex')
legend('$\Upsilon_{B}$ $(N)$','$\Upsilon_{D}$ $(N)$','$\Upsilon_{F}$ $(N)$','Interpreter','latex','FontSize', font_size, 'Location', 'northwest')
axis([0 25206.93764 -20*10^-10 20*10^-10])
    ax = gca;
    ax.FontSize = font_size;
    ax = gca; ax.FontSize = 25;

figure()
% grid on;
   % Set figure dimensions (width x height in pixels)
    fig_width = 600;
    fig_height = 500;
    set(gcf, 'Position', [100, 100, fig_width, fig_height]);
    % Set font size (customizable)
    font_size = 16.5; % Change this value as needed
plot( someName(:,1), someName(:,2), 'b', someName(:,1), someName(:,3), 'k',someName(:,1), someName(:,4), 'r' )
xlabel('t (s)'); ylabel('States', 'Interpreter','latex');
hold on
% grid on;
plot(ts,SV(:,1),'--') 
legend('$q_{7}$ (rad)','$q_{8}$ (m)','$q_{9}$ (m)','$q_{8_{d}}$ (m)','Interpreter','latex','FontSize', font_size, 'Location', 'northwest')
axis([0 25206.93764 -5*10^-6 60*10^-6])
    ax = gca;
    ax.FontSize = font_size;
    ax = gca; ax.FontSize = 25;

% figure(3)
% grid on;
% plot( someName(:,1), someName(:,2), 'b', someName(:,1), someName(:,3), 'k',someName(:,1), someName(:,4), 'r' )
% xlabel('t (s)'); ylabel('States', 'Interpreter','latex');
% legend('$q_{7}$ (rad)','$q_{8}$ (m)','$q_{9}$ (m)','Interpreter','latex')


%ax = gca;
%ax.FontSize = 30; 
%exportgraphics(ax,'CellSimplifedModV3_Version3_OptimalControl.pdf','BackgroundColor','none')