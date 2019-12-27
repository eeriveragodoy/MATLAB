%% Edgardo Rivera Godoy
%  912753270
%  Project 9

function [] = prob2(Re)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Grid
dx = 0.01; % x spacing
dy = 0.01; % y spacing
x = -2*dx:dx:4; % x range
y = 0:dy:1; % y range

sx = size(x); % number of x points
sy = size(y); % number of y points

LE = find(x == 0); % leading edge location
TE = find(x == 1); % trailing edge location

[Y,X] = meshgrid(y,x); % creates X Y meshgrid
u = ones(size(X)); % initial velocity guess
v = zeros(size(X)); % initial velocity guess

tau = zeros(sx); % shear stress on surface of plate

% First guess for tridiagonal matrix coefficients: a, b, c, and d
a = zeros(sy);
b = a;
c = a;
d = a;
b(end) = 1;
d(end) = 1;

j = 2:length(y)-1;
for i = 2:length(x)   %i <-- 2 to imax
    u(i,:) = u(i-1,:);
    v(i,:) = v(i-1,:);
    
    %Compute tridiagonal matrix coefficients: a, b, c, and d
    a(j) = -1/(Re*dy^2) - v(i,j)/(2*dy);
    b(j) =  2/(Re*dy^2) + u(i,j)/dx;
    c(j) = -1/(Re*dy^2) + v(i,j)/(2*dy);
    d(j) = u(i-1,j).*u(i,j)./dx;
    
    if i <= TE % before trailing edge
        b(1) = 1;
        c(1) = 0;
        d(1) = 0;
    else  % after trailing edge
        b(1) = u(i,1)/dx + 2/(Re*dy^2);
        c(1) = -2/(Re*dy^2);
        d(1) = (u(i,1) * u(i-1,1))/dx;
    end
    % Define
    iter = 1;
    res = 1e8;
    tol = 1e-6; % max allowed error
    
    while res > tol
        u_star = thomas3(a,b,c,d); % tridiagonal solver
        u(i,:) = u_star; % use newest solution
        % obtain new v
        if  i > LE
            v(i,j) = v(i,j-1) - dy/(2*dx)*(u(i,j) - u(i-1,j) + u(i,j-1) - u(i-1,j-1));
        end
        % recalculate coefficients
        a(j) = -1/(Re*dy^2) - v(i,j)/(2*dy);
        b(j) =  2/(Re*dy^2) + u(i,j)/dx;
        c(j) = -1/(Re*dy^2) + v(i,j)/(2*dy);
        d(j) = u(i-1,j).*u(i,j)/dx;
        if i <= TE
            b(1) = 1;
            c(1) = 0;
            d(1) = 0;
        else 
            b(1) = u(i,1)/dx + 2/(Re*dy^2);
            c(1) = -2/(Re*dy^2);
            d(1) = (u(i,1) * u(i-1,1))/dx;
        end
        % calculate residuals
        res = max(abs(a(j).*u(i,j-1) + b(j).*u(i,j) + c(j) .* u(i,j+1) - d(j)));
        res = max(res, abs(b(1)*u(i,1) + c(1)*u(i,2) - d(1)));
        iter = iter + 1;
    end
    u_y = (-u(i,3) + 4*u(i,2) - 3*u(i,1))/(2*dx); % du/dy
    tau(i) = 2.*u_y./Re; % shear stress at x(i)
end

% countour of u plot
figure()
contourf(X,Y,u,30,'linecolor','none')
colormap jet
colorbar;
title('Contour of $u$', 'interpreter', 'latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

xs = [0, 0.25, 0.5, 0.75, 1.0, 1+dx,2,3,4];
index = zeros(1,length(xs));
for i = 1:length(xs)
    index(i) = find(x == xs(i));
end


figure();
hold on
% Velocity distribution
s = quiver(X(index,:),Y(index,:),u(index,:),zeros(length(index),length(y)), 'linewidth', 1);
s.MaxHeadSize = 0.003; 
% Lamninar BL thickness
delta = 5.*x(LE:TE)./sqrt(Re.*x(LE:TE)); 
plot(x(LE:TE), delta', 'r', 'linewidth', 2) 
hold on
% Flat Plate
plot(x(LE:TE), zeros(size(x(LE:TE))), 'k', 'linewidth', 3);
legend('U(y)', 'Laminar BL Thickness', 'Plate');
title(['BL on a Flat Plate at Re = ' num2str(Re)])
xlabel('x')
ylabel('y')
hold off

% Shear stress on surface plot
figure()
plot(x,tau,'k','linewidth',2)
xlim([-dx,1])
title('Shear Stress on the surface of the plate','interpreter','latex')
xlabel('x','interpreter','latex')
ylabel('$\tau$','interpreter','latex')
grid on
hold off

end

