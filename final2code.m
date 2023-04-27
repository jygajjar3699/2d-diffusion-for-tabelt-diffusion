% Define parameters
D = 2e-9; % Diffusion coefficient of the dissolved substance
L = 0.01; % Length of the domain
N = 100; % Number of grid points
dx = L/N; % Grid spacing
C0 = 10; % Concentration of the dissolved substance in the tablet
h = 0.001; % Height of the tablet
w = 0.005; % Width of the tablet
x0 = 0.005; % x-coordinate of the center of the tablet
y0 = 0.01; % y-coordinate of the center of the tablet
R = 0.001; % Radius of the tablet (assumed circular)
v = 0.002; % Velocity of the tablet
t_end = 100; % End time of simulation
dt = 0.1*dx^2/D; % Time step
rho_fluid = 1000; % Density of the liquid
rho_tablet = 2000; % Density of the tablet
% Create grid of x- and y-coordinates
x = linspace(0, L, N);
y = linspace(0, L, N);
[X, Y] = meshgrid(x, y);
% Initialize concentration and density matrices
C = zeros(N, N);
rho = rho_fluid * ones(N, N);
% Set initial concentration and density
C((Y - y0).^2 + (X - x0).^2 <= R^2) = C0;
% rho(idx) = rho_tablet;
% Create a new VideoWriter object
writerObj = VideoWriter('finaloutput.mp4', 'MPEG-4');
writerObj.FrameRate = 30; % set the frame rate of the video
open(writerObj); % open the video writer object
% Create figure for animation
fig = figure;
set(fig, 'Position', [100, 100, 800, 500]);
ax = gca;
ax.NextPlot = 'replace Children';
ax.XLim = [0 L];
ax.YLim = [0 L];
ax.ZLim = [0 C0];
ax.XLabel.String = 'x (m)';
ax.YLabel.String = 'y (m)';
ax.ZLabel.String = 'Concentration';
% Simulate diffusion and advection
frames=[];
for t = dt:dt:t_end
 % Update tablet position
 y0 = y0 + v*dt;
 % Compute concentration and density at next time step
 C_new = C + D*dt/dx^2*(circshift(C, [0 -1]) + circshift(C, [0 1]) + circshift(C, [-1 0]) + circshift(C, [1 0]) - 4*C);
 rho_new = rho_fluid + (rho_tablet - rho_fluid) * (C_new > 0);
 rho_new = rho_new + v*dt/dx*(circshift(rho_new, [0 -1]) - rho_new) + D*dt/dx^2*(circshift(rho_new, [0 -1]) + circshift(rho_new, [0 1]) + circshift(rho_new, [-10]) + circshift(rho_new, [1 0]) - 4*rho_new);
 % Apply boundary conditions for concentration (no flux)
 C_new(1,:) = C_new(2,:);
 C_new(end,:) = C_new(end-1,:);
 C_new(:,1) = C_new(:,2);
 C_new(:,end) = C_new(:,end-1);

 % Apply boundary conditions for density (no flux at top and bottom, free surface at sides)
 rho_new(1,:) = rho_new(2,:);
 rho_new(end,:) = rho_new(end-1,:);
 rho_new(:,1) = rho_new(:,2);
 rho_new(:,end) = rho_new(:,end-1);


 % Update concentration and density matrices
 C = C_new;
 rho = rho_new;
 % Create frame for current time step
 frames(end+1) = surf(ax, X, Y, C,rho);
 title(ax, sprintf('Time: %.2f s', t));

 % Pause to control animation speed
 pause(0.01);
 % Write the current frame to the video
 frame = getframe(gcf); % capture the current figure as a frame
 writeVideo(writerObj, frame); % write the frame to the video
end
close(writerObj); % close the video writer object