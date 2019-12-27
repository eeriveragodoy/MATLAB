%% Edgardo Rivera Godoy
%  912753270
%  Project 9
%% Problem I 2D compressible inviscid flow over a thin airfoil 
% The supersonic flow over a thin airfoil is governed by the following:
% (1 - M_inf^2)*u_xx + u_yy = 0, M_inf = 0.8, 1.4, and 1.8; u_y = v_x
% Where u is the normalized x component of the flow and v is the normalized y component of the flow. 
% Use the following initial and boundary conditions (note * only applies to subsonic case):
% (-1, 0) to (-1, 4): u = 1
% (-1, 4) to (2, 4): u = 1
% (2, 4) to (2, 0): u = 1*
% (-1, 0) to (0, 0): v = 0 
% (0, 0) to (1, 0): v = dyB2dx
% (1, 0) to (2, 0): v = 0 
% dx = dy = 0.02
% Plot the contours of the pressure coefficient for each case. Cp = ?2(u ? 1) 

% Grid
dx = 0.02; % x spacing
dy = 0.02; % y spacing
x = -1:dx:2; % x range
y = 0:dy:4; % y range
[X, Y] = meshgrid(x, y); % X, y meshgrids
M = 0.8; % subsonic mach number

% Airfoil
tau = 1/10;
yB = [zeros(size(x(x < 0))), tau*sin(pi*x(x>=0 & x <=1)), zeros(size(x(x > 1)))];
dyB2dx = [0, yB(3:end) - 2.*yB(2:end-1) + yB(1:end-2), 0]./dx^2;

% Iterations
% tridiagonal system coefficients a b c
a = [0; 1/dy^2 * ones(length(y)-2,1); 0];
c = [2/dy^2; 1/dy^2 * ones(length(y)-2,1); 0];
b = [-2.*(1/dy^2 + (1-M.^2)/dx^2) * ones(length(y)-1,1); 1];

u = ones(size(X)); % orignal guess for solution
err = 1e-5; % min allowed error
w = 1.9;
res = 1000;

while res(end) > err % loops while latest error is too high
    res(end+1) = 0;
    for i = 2:length(x)-1 % Computes tridiagonal matrix coefficient: d
       d = -(1-M.^2).*(u(:,i+1) + u(:,i-1))./dx^2;
       d(1) = d(1) + 2*dyB2dx(i)/dy;
       d(end) = u(end,i);
       
       % Gauss Elimination for Tridiagonal Systems
       u_star = thomas3(a, b, c, d);
       u(:,i) = w*u_star + (1-w).*u(:,i);

       % calculate residuals
       resIter = abs([  b(1).*u(1,i) + c(1).*u(2,i) - d(1); a(2:end-1).*u(1:end-2,i) + ...
                        b(2:end-1).*u(2:end-1,i) + c(2:end-1).*u(3:end,i) - d(2:end-1)]);
       if max(resIter) > res(end)
           res(end) = max(resIter); % max residual is the current iterations residual
       end
    end
end
% pressure coefficient
Cp = 2*(1 - u);

% Contour plots of the pressure coefficient
figure();
contourf(X, Y, u, 20);
colormap winter
h = colorbar;
ylabel(h, 'u')
title('Contour of $u$, $M_\infty = 0.8$', 'interpreter', 'latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

% Plot Pressure Coefficient
figure();
contourf(X, Y, Cp, 20);
colormap jet
h = colorbar;
ylabel(h, 'c_p')
title('Contours of Pressure Coefficient, $M_\infty = 0.8$', 'interpreter', 'latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

% Residual
figure();
lres = 1:length(res);
plot(lres, res,'linewidth',2);
set(gca, 'yscale', 'log')
grid on
title('Residual over Iterations, $M_\infty = 0.8$', 'interpreter', 'latex');
xlabel('$Iteration$','interpreter','latex');
ylabel('$Residual$','interpreter','latex');

%% M = 1.4
% Grid
dx = 0.02;
dy = 0.02;
x = -1:dx:2;
y = 0:dy:4;
[X, Y] = meshgrid(x, y);

lx = size(X);
M = 1.4;
beta = (1 - M^2);

% Airfoil
tau = 1/10;
yB = [zeros(size(x(x < 0))), tau*sin(pi*x(x>=0 & x <=1)), zeros(size(x(x > 1)))];
dyB2dx = [0, yB(3:end) - 2.*yB(2:end-1) + yB(1:end-2), 0]./dx^2;

% Iterations
a = [0; 1/dy^2 * ones(length(y)-2,1); 0];
c = [2/dy^2; 1/dy^2 * ones(length(y)-2,1); 0];
b = [beta/dx^2 - 2/dy^2; ( beta/dx^2 - 2/dy^2 ) * ones(length(y)-2,1); 1];

u = ones(length(x),length(y));

a = ones(1, length(y))*1/dy^2; 
b = ones(1, length(y))*(beta/dx^2 - 2/dy^2); 
c = ones(1, length(y))*1/dy^2; 
d  = zeros(1, length(y)); 

c(1) = 2/dy^2;
b(end) = 1;
a(1) = 0; 
a(end) = 0; 
d(end) = 1; 

for i = 3:length(x)-1
    for j = 2:length(y) - 1
        d(j) = (-u(i-2,j) + 2*u(i-1,j))*beta / dx^2;
    end
    d(1) = (-u(i-2,1) + 2*u(i-1,1))*beta / dx^2 + 2/dy*dyB2dx(i);
    d(end) = u(i, end);
    utemp = thomas3(a, b, c, d);
    u(i, :) =utemp; 
end
    
Cp = 2*(1 - u);

% Plot Contour
figure();
contourf(X', Y', u, 20);
colormap winter
h = colorbar;
ylabel(h, 'u')
title('Contour of $u$, $M_\infty = 1.4$', 'interpreter', 'latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

% Plot Pressure Coefficient
figure();
% plot(x, Cp(:,1),'linewidth',2);
% set(gca,'Ydir','reverse')
% grid on
contourf(X', Y', Cp, 20);
colormap winter
h = colorbar;
ylabel(h, 'c_p')
title('Contours of Pressure Coefficient, $M_\infty = 1.4$', 'interpreter', 'latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

%% M = 1.8
% Grid
dx = 0.02;
dy = 0.02;
x = -1:dx:2;
y = 0:dy:4;
[X, Y] = meshgrid(x, y);

lx = size(X);
M = 1.8;
beta = (1 - M^2);

% Airfoil
tau = 1/10;
yB = [zeros(size(x(x < 0))), tau*sin(pi*x(x>=0 & x <=1)), zeros(size(x(x > 1)))];
dyB2dx = [0, yB(3:end) - 2.*yB(2:end-1) + yB(1:end-2), 0]./dx^2;

% Iterations
a = [0; 1/dy^2 * ones(length(y)-2,1); 0];
c = [2/dy^2; 1/dy^2 * ones(length(y)-2,1); 0];
b = [beta/dx^2 - 2/dy^2; ( beta/dx^2 - 2/dy^2 ) * ones(length(y)-2,1); 1];

u = ones(length(x),length(y));

a = ones(1, length(y))*1/dy^2; 
b = ones(1, length(y))*(beta/dx^2 - 2/dy^2); 
c = ones(1, length(y))*1/dy^2; 
d  = zeros(1, length(y)); 

c(1) = 2/dy^2;
b(end) = 1;
a(1) = 0; 
a(end) = 0; 
d(end) = 1; 

for i = 3:length(x)-1
    for j = 2:length(y) - 1
        d(j) = (-u(i-2,j) + 2*u(i-1,j))*beta / dx^2;
    end
    d(1) = (-u(i-2,1) + 2*u(i-1,1))*beta / dx^2 + 2/dy*dyB2dx(i);
    d(end) = u(i, end);
    utemp = thomas3(a, b, c, d);
    u(i, :) =utemp;

end
    
Cp = 2*(1 - u);

% Plot Contour
figure();
contourf(X', Y', u, 20);
colormap winter
h = colorbar;
ylabel(h, 'u')
title('Contour of $u$, $M_\infty = 1.8$', 'interpreter', 'latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

% Plot Pressure Coefficient
figure();
% plot(x, Cp(:,1),'linewidth',2);
% set(gca,'Ydir','reverse')
% grid on
contourf(X', Y', Cp, 20);
colormap winter
h = colorbar;
ylabel(h, 'c_p')
title('Pressure Coefficient, $M_\infty = 1.8$', 'interpreter', 'latex');
xlabel('$x$','interpreter','latex');
ylabel('$C_p$','interpreter','latex');


%% Problem 2 Incompressible viscous boundary layer flow over a flat plate 
% The viscous flow over a flat plate is governed by the boundary layer equations: 
% u_x + v_y = 0; uu_x + vu_y = u_yy/Re_L
% Use the following initial and boundary conditions (note * only applies to subsonic case):
% (-2dx, 0) to (-2dx, 1): u = 1
% (-2dx, 1) to (4, 1): u = 1
% (4, 1) to (4, 0):
% (-2dx, 0) to (0, 0): v = 0, u_y = 0
% (0, 0) to (1, 0): u = v = 0
% (1, 0) to (4, 0): v = 0, u_y =0
% dx = dy = 0.01
% Find u and v for Re_L ? {100,1000,10000} Plot the velocity contours. 
% Plot u(x,y) vs y (u on the x-axis) at x ? {0.0,0.25,0.5,0.75,1.0,1 + ?x,2.0,3.0,4.0} 
% On the same plot, plot the boundary layer thickness over the plate (0 ? x ? 1). 
% The boundary layer thickness is approximated by:
% delta/x = 5/sqrt(Re_x), Re_x = rho*U_inf*x/mu
% Plot the shear stress on the surface of the plate using a 2nd order accurate finite difference scheme 
 
prob2(100)      % Re = 100

prob2(1000)     % Re = 1000

prob2(10000)    % Re = 10000