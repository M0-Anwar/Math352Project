function project()
    simple_fluid_gui();
end

function simple_fluid_gui()
    % Create the main figure
    f = figure('Name', 'Simple Fluid Simulator', ...
               'Position', [100, 100, 1000, 700], ...
               'Color', [0.1 0.1 0.1]);  % Dark background

    % --- UI Controls (Inputs) ---
    uicontrol('Style', 'text', 'Position', [40 360 100 20], 'String', 'Viscosity:', ...
              'BackgroundColor', [0.1 0.1 0.1], 'ForegroundColor', [1 1 1]);
    viscosity_edit = uicontrol('Style', 'edit', 'Position', [150 360 120 25], 'String', '1', ...
              'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1]);

    uicontrol('Style', 'text', 'Position', [40 320 100 20], 'String', 'Pressure K:', ...
              'BackgroundColor', [0.1 0.1 0.1], 'ForegroundColor', [1 1 1]);
    K_edit = uicontrol('Style', 'edit', 'Position', [150 320 120 25], 'String', '500', ...
              'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1]);

    uicontrol('Style', 'text', 'Position', [40 280 100 20], 'String', 'Blob X center:', ...
              'BackgroundColor', [0.1 0.1 0.1], 'ForegroundColor', [1 1 1]);
    xcenter_edit = uicontrol('Style', 'edit', 'Position', [150 280 120 25], 'String', '4', ...
              'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1]);

    uicontrol('Style', 'text', 'Position', [40 220 100 20], 'String', 'Blob Y center:', ...
              'BackgroundColor', [0.1 0.1 0.1], 'ForegroundColor', [1 1 1]);
    ycenter_edit = uicontrol('Style', 'edit', 'Position', [150 220 120 25], 'String', '128', ...
              'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1]);

    uicontrol('Style', 'text', 'Position', [40 180 100 20], 'String', 'Blob Radius:', ...
              'BackgroundColor', [0.1 0.1 0.1], 'ForegroundColor', [1 1 1]);
    radius_edit = uicontrol('Style', 'edit', 'Position', [150 180 120 25], 'String', '5', ...
              'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1]);

    start_button = uicontrol('Style', 'pushbutton', 'String', 'Start Simulation', ...
                            'Position', [40 120 230 40], ...
                            'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1], ...
                            'Callback', @start_simulation);

    % --- Axes for plotting ---
    ax = axes('Units', 'pixels', 'Position', [300, 50, 600, 650], ...
              'XTick', [], 'YTick', [], 'Box', 'on', 'Color', [0 0 0]);

    running = false; % Control flag

    % --- Callback function ---
    function start_simulation(~,~)
        if running
            return; % Don't start multiple times
        end
        running = true;
        simulate_fluid(ax, viscosity_edit, K_edit, xcenter_edit, ycenter_edit, radius_edit);
    end
end

function simulate_fluid(ax, viscosity_edit, K_edit, xcenter_edit, ycenter_edit, radius_edit)
    % Grid settings
    Nx = 256; Ny = 256;
    dx = 1.0; dy = 1.0;
   dt = 0.01;
  % max_vel = max(max(sqrt(u.^2 + v.^2)));
%dt = 0.25 * dx / (max_vel + 1e-5);  % Add small number to avoid division by zero


    % Initialize fields
    %u = zeros(Nx, Ny);
    %v = zeros(Nx, Ny);
    %density = zeros(Nx, Ny);
    % Initialize fields
  u = 0.0 *abs((zeros(Nx))); % Initial velocity to the right (positive x direction)
  v = 0.0*abs(sin(zeros(Ny)));      % No vertical movement
  density = sin(zeros(Nx, Ny));


    [X, Y] = meshgrid(1:Nx, 1:Ny);
    cx = Nx/2; cy = Ny/2;  % Center of vortex

% Rotational velocity field
u = -(Y - cy) ./ sqrt((X - cx).^2 + (Y - cy).^2 + 1e-5);
v =  (X - cx) ./ sqrt((X - cx).^2 + (Y - cy).^2 + 1e-5);

strength = 5.0;
u = strength * u;
v = strength * v;


    % Setup plotting
    axes(ax);
    colormap('jet');

    theta =0 ;
    % Simulation Loop
    for t = 1:10000
        % --- Read Parameters Live ---
        viscosity = str2double(get(viscosity_edit, 'String'));

        K = str2double(get(K_edit, 'String'));


        x_center = str2double(get(xcenter_edit, 'String'));
        y_center = str2double(get(ycenter_edit, 'String'));
        radius = str2double(get(radius_edit, 'String'));







        %radius_motion = 50;
         %theta = theta + 0.01;
          %x_center = 128 +  10 * cos(theta)  ;
          %y_center = + (sin(theta) + cos(theta)) ;


        % === 1. Re-create density blob if needed ===
        blob = ((X-x_center).^2)/0.5 + ((Y-y_center).^2)/9< radius^2;
        density(blob) = max(density(blob),0.8);

        % === 2. Compute Derivatives ===
        grad_rho_x = (circshift(density, [-1,0]) - circshift(density, [1,0])) / (2*dx);
        grad_rho_y = (circshift(density, [0,-1]) - circshift(density, [0,1])) / (2*dy);
        div_u = (circshift(u, [-1,0]) - circshift(u, [1,0])) / (2*dx) + ...
                (circshift(v, [0,-1]) - circshift(v, [0,1])) / (2*dy);

        % === 3. Update Density ===
        density = density + dt * (-u .* grad_rho_x - v .* grad_rho_y - density .* div_u);
        %density = max(density, 0);
        density = min(max(density, 0), 1.0);
        % Zero density at boundaries
density(1,:)   = 0;   % Top row
density(end,:) = 0;   % Bottom row
density(:,1)   = 0;   % Left column
density(:,end) = 0;   % Right column



        % === 4. Semi-Lagrangian Advection for velocity ===
        [Xq, Yq] = meshgrid(1:Nx, 1:Ny);
        X_back = Xq - dt * u;
        Y_back = Yq - dt * v;
        X_back = min(max(X_back, 1), Nx);
        Y_back = min(max(Y_back, 1), Ny);
        u = interp2(Xq, Yq, u, X_back, Y_back, 'linear', 0);
        v = interp2(Xq, Yq, v, X_back, Y_back, 'linear', 0);

        % === 5. Forces ===
        gradP_x = K * grad_rho_x;
        gradP_y = K * grad_rho_y;
        Lap_u = (circshift(u, [-1,0]) + circshift(u, [1,0]) + circshift(u, [0,-1]) + circshift(u, [0,1]) - 4*u) / (dx*dx);
        Lap_v = (circshift(v, [-1,0]) + circshift(v, [1,0]) + circshift(v, [0,-1]) + circshift(v, [0,1]) - 4*v) / (dy*dy);
        u = u + dt * (-gradP_x + viscosity * Lap_u);
        v = v + dt * (-gradP_y + viscosity * Lap_v);
        u = u + 0.01 * (2*rand(Nx,Ny) - 1);
        v = v + 0.01 * (2*rand(Nx,Ny) - 1);

        % === 6. Boundary Conditions ===
        u([1 end],:) = 0; u(:,[1 end]) = 0;
        v([1 end],:) = 0; v(:,[1 end]) = 0;

        % === 7. Visualization ===
        if mod(t,10) == 0
            imagesc(density',[0 1]);
            axis equal tight off;
            title(['Time step: ', num2str(t)]);
            colorbar;
            drawnow;

        end

        pause(0); % Pause for GUI responsiveness
    end % loop
end % function
