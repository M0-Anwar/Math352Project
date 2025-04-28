function project()
    simple_fluid_gui();
end

function simple_fluid_gui()
    % Create the main figure
    f = figure('Name', 'Simple Fluid Simulator', ...
               'Position', [100, 100, 1000, 600], ...
               'Color', [0.1 0.1 0.1]);  % Dark background

    % --- UI Controls (Inputs) ---
    uicontrol('Style', 'text', 'Position', [40 360 100 20], 'String', 'Viscosity:', ...
              'BackgroundColor', [0.1 0.1 0.1], 'ForegroundColor', [1 1 1]);
    viscosity_edit = uicontrol('Style', 'edit', 'Position', [150 360 120 25], 'String', '0.001', ...
              'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1]);

    uicontrol('Style', 'text', 'Position', [40 320 100 20], 'String', 'Pressure K:', ...
              'BackgroundColor', [0.1 0.1 0.1], 'ForegroundColor', [1 1 1]);
    K_edit = uicontrol('Style', 'edit', 'Position', [150 320 120 25], 'String', '0.5', ...
              'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1]);

    uicontrol('Style', 'text', 'Position', [40 280 100 20], 'String', 'Blob X center:', ...
              'BackgroundColor', [0.1 0.1 0.1], 'ForegroundColor', [1 1 1]);
    xcenter_edit = uicontrol('Style', 'edit', 'Position', [150 280 120 25], 'String', '32', ...
              'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1]);

    uicontrol('Style', 'text', 'Position', [40 220 100 20], 'String', 'Blob Y center:', ...
              'BackgroundColor', [0.1 0.1 0.1], 'ForegroundColor', [1 1 1]);
    ycenter_edit = uicontrol('Style', 'edit', 'Position', [150 220 120 25], 'String', '32', ...
              'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1]);

    uicontrol('Style', 'text', 'Position', [40 180 100 20], 'String', 'Blob Radius:', ...
              'BackgroundColor', [0.1 0.1 0.1], 'ForegroundColor', [1 1 1]);
    radius_edit = uicontrol('Style', 'edit', 'Position', [150 180 120 25], 'String', '10', ...
              'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1]);

    start_button = uicontrol('Style', 'pushbutton', 'String', 'Start Simulation', ...
                            'Position', [40 120 230 40], ...
                            'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1], ...
                            'Callback', @start_simulation);

    % --- Axes for plotting ---
    ax = axes('Units', 'pixels', 'Position', [350, 50, 500, 500], ...
              'XTick', [], 'YTick', [], 'Box', 'on', 'Color', [0 0 0]);

    % --- Callback function ---
    function start_simulation(~,~)
        viscosity = str2double(get(viscosity_edit, 'String'));
        K = str2double(get(K_edit, 'String'));
        x_center = str2double(get(xcenter_edit, 'String'));
        y_center = str2double(get(ycenter_edit, 'String'));
        radius = str2double(get(radius_edit, 'String'));

        simulate_fluid(viscosity, K, x_center, y_center, radius, ax);
    end
end


function simulate_fluid(viscosity, K, x_center, y_center, radius, ax)
    % Grid settings
    Nx = 64; Ny = 64;
    dx = 1.0; dy = 1.0;
    dt = 0.1;

    % Initialize fields
    u = zeros(Nx, Ny);
    v = zeros(Nx, Ny);
    density = zeros(Nx, Ny);

    [X, Y] = meshgrid(1:Nx, 1:Ny);
    density( (X-x_center).^2 + (Y-y_center).^2 < radius^2 ) = 1.0;

    % Setup plotting
    axes(ax);
    colormap('jet');

    % Simulation Loop
    for t = 1:500
        % === 1. Compute Derivatives ===
        grad_rho_x = (circshift(density, [-1,0]) - circshift(density, [1,0])) / (2*dx);
        grad_rho_y = (circshift(density, [0,-1]) - circshift(density, [0,1])) / (2*dy);
        div_u = (circshift(u, [-1,0]) - circshift(u, [1,0])) / (2*dx) + ...
                (circshift(v, [0,-1]) - circshift(v, [0,1])) / (2*dy);

        % === 2. Update Density ===
        density = density + dt * (-u .* grad_rho_x - v .* grad_rho_y - density .* div_u);
        density = max(density, 0);

        % === 3. Semi-Lagrangian Advection for velocity ===
        [Xq, Yq] = meshgrid(1:Nx, 1:Ny);
        X_back = Xq - dt * u;
        Y_back = Yq - dt * v;
        X_back = min(max(X_back, 1), Nx);
        Y_back = min(max(Y_back, 1), Ny);
        u = interp2(Xq, Yq, u, X_back, Y_back, 'linear', 0);
        v = interp2(Xq, Yq, v, X_back, Y_back, 'linear', 0);

        % === 4. Forces ===
        gradP_x = K * grad_rho_x;
        gradP_y = K * grad_rho_y;
        Lap_u = (circshift(u, [-1,0]) + circshift(u, [1,0]) + circshift(u, [0,-1]) + circshift(u, [0,1]) - 4*u) / (dx*dx);
        Lap_v = (circshift(v, [-1,0]) + circshift(v, [1,0]) + circshift(v, [0,-1]) + circshift(v, [0,1]) - 4*v) / (dy*dy);
        u = u + dt * (-gradP_x + viscosity * Lap_u);
        v = v + dt * (-gradP_y + viscosity * Lap_v);

        % === 5. Boundary Conditions ===
        u([1 end],:) = 0; u(:,[1 end]) = 0;
        v([1 end],:) = 0; v(:,[1 end]) = 0;

        % === 6. Visualization ===
        if mod(t,10) == 0
            imagesc(density');
            axis equal tight;
            title(['Time step: ', num2str(t)]);
            colorbar;
            drawnow;
        end % if
    end % loop
end % function

