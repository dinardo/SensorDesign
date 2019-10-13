%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the velocity field taking into account      %
% the magnetic field (coming out of the surface/screen) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relations between Force field and velocity field
% mu * (Ex + vy*B) = vx --> vx = mu / (1 + (mu*R*B)^2) * (Ex + mu*R*B*Ey)
% --> vx = mu / [(1+(mu*R*B)^2) * (1+mu*E/vs)] * (Ex + mu*R*B*Ey)
% mu * (Ey - vx*B) = vy --> vy = mu / (1 + (mu*R*B)^2) * (Ey - mu*R*B*Ex)
% --> vy = mu / [(1+(mu*R*B)^2) * (1+mu*E/vs)] * (Ey - mu*R*B*Ex)

% Potential = Solution of the Poisson equation
% Step      = Unit step of the lattice on which the field is computed [um]
% Bulk      = Bulk thickness [um]
% BField    = Magnetic field (orthogonal+outgoing from the 2D geometry) [T]
% Pitch     = Strip pitch [um]

% mu_e      = Electron mobility [um^2/(V*ns)]
% RH_e      = Relative Hall electron mobility
% vs_e      = Saturation velocity of the electrons [um/ns]
% beta_e    = Exponent for the electric field dependence of mobility

% mu_h      = Hole mobility in [um^2/(V*ns)]
% RH_h      = Relative Hall hole mobility
% vs_h      = Saturation velocity of the holes [um/ns]
% beta_h    = Exponent for the electric field dependence of mobility

% ItFigIn   = Figure iterator input

function [VFieldx_e, VFieldy_e, VFieldx_h, VFieldy_h, x, y, ItFigOut] =...
    VelocityField(Potential,Step,Bulk,BField,Pitch,...
    mu_e,RH_e,vs_e,beta_e,mu_h,RH_h,vs_h,beta_h,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReSample   = 2;   % Used in order to make nice plots [um]
ContLevel  = 40;  % Contour plot levels
MagnVector = 1.5; % Vector field magnification

x = -Pitch:Step:Pitch; % Bound along x-coordinate
y = 0:Step:Bulk * 3/2; % Bound along y-coordinate

VFieldx_e = zeros(length(y), length(x));
VFieldy_e = zeros(length(y), length(x));
VFieldx_h = zeros(length(y), length(x));
VFieldy_h = zeros(length(y), length(x));


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
[mshx,mshy] = meshgrid(x,y);
queryPoints = [mshx(:),mshy(:)]';

[Etx, Ety] = evaluateGradient(Potential,queryPoints);

Etx = reshape(Etx,size(mshx));
Ety = reshape(Ety,size(mshy));

interp = interpolateSolution(Potential,queryPoints);
interp = reshape(interp,size(mshx));


% Velocity field calculation for electrons
fprintf('@@@ I''m calculating the Velocity-Field @@@\n');
for i = 1:length(x)
    for j = 1:length(y)
        
        E = sqrt(Etx(j,i)^2 + Ety(j,i)^2);

        if ~isnan(Etx(j,i)) && ~isnan(Ety(j,i))
            q = -1; % Electron charge sign
            VFieldx_e(j,i) = mu_e / ((1 + (mu_e*RH_e*BField*1e-3)^2) *...
                (1 + (mu_e*E/vs_e)^beta_e)^(1/beta_e)) *...
                (q*Etx(j,i) + mu_e*RH_e*BField*1e-3*Ety(j,i));
            
            VFieldy_e(j,i) = mu_e / ((1 + (mu_e*RH_e*BField*1e-3)^2) *...
                (1 + (mu_e*E/vs_e)^beta_e)^(1/beta_e)) *...
                (q*Ety(j,i) - mu_e*RH_e*BField*1e-3*Etx(j,i));

            
            q = +1; % Hole charge sign
            VFieldx_h(j,i) = mu_h / ((1 + (mu_h*RH_h*BField*1e-3)^2) *...
                (1 + (mu_h*E/vs_h)^beta_h)^(1/beta_h)) *...
                (q*Etx(j,i) + mu_h*RH_h*BField*1e-3*Ety(j,i));

            VFieldy_h(j,i) = mu_h / ((1 + (mu_h*RH_h*BField*1e-3)^2) *...
                (1 + (mu_h*E/vs_h)^beta_h)^(1/beta_h)) *...
                (q*Ety(j,i) - mu_h*RH_h*BField*1e-3*Etx(j,i));
        else
            VFieldx_e(j,i) = 0;
            VFieldy_e(j,i) = 0;

            VFieldx_h(j,i) = 0;
            VFieldy_h(j,i) = 0; 
        end
    end
end

% Lorentz angle in the middle of the strip on the backplane side
Lore = abs(atan(VFieldx_e(1,int32(length(x)/2)) /...
    VFieldy_e(1,int32(length(x)/2))) * 180/pi);
fprintf('Lorentz angle for electrons: %0.1f [degree]\n',Lore);
Lore = abs(atan(VFieldx_h(1,int32(length(x)/2)) /...
    VFieldy_h(1,int32(length(x)/2))) * 180/pi);
fprintf('Lorentz angle for holes: %0.1f [degree]\n',Lore);


%%%%%%%%%
% Plots %
%%%%%%%%%
VFieldx_e_ReSample = VFieldx_e(1:ReSample:length(y), 1:ReSample:length(x));
VFieldy_e_ReSample = VFieldy_e(1:ReSample:length(y), 1:ReSample:length(x));
VFieldx_h_ReSample = VFieldx_h(1:ReSample:length(y), 1:ReSample:length(x));
VFieldy_h_ReSample = VFieldy_h(1:ReSample:length(y), 1:ReSample:length(x));

xp = x(1:ReSample:length(x));
yp = y(1:ReSample:length(y));

[xx,yy] = meshgrid(x,y);

figure(ItFigIn);
colormap jet;
subplot(2,2,1);
surf(xx,yy,VFieldx_e,'EdgeColor','none','FaceColor','interp');
title('X component of the Velocity-Field for e');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('Velocity Field [\mum / ns]');
subplot(2,2,2);
surf(xx,yy,VFieldy_e,'EdgeColor','none','FaceColor','interp');
title('Z component of the Velocity-Field for e');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('Velocity Field [\mum / ns]');
subplot(2,2,3);
surf(xx,yy,VFieldx_h,'EdgeColor','none','FaceColor','interp');
title('X component of the Velocity-Field for h');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('Velocity Field [\mum / ns]');
subplot(2,2,4);
surf(xx,yy,VFieldy_h,'EdgeColor','none','FaceColor','interp');
title('Z component of the Velocity-Field for h');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('Velocity Field [\mum / ns]');

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
colormap jet;
subplot(1,2,1);
contour(x,y,interp,ContLevel);
hold on
quiver(xp,yp,VFieldx_e_ReSample,VFieldy_e_ReSample,MagnVector);
colormap jet;
hold off
title('Velocity-Field for electrons');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('Velocity Field [\mum / ns]');
subplot(1,2,2);
contour(x,y,interp,ContLevel);
hold on
quiver(xp,yp,VFieldx_h_ReSample,VFieldy_h_ReSample,MagnVector);
colormap jet;
hold off
title('Velocity-Field for holes');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('Velocity Field [\mum / ns]');

ItFigOut = ItFigIn + 1;
fprintf('CPU time --> %.2f [min]\n\n',(cputime-TStart)/60);
end
