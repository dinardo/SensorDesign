%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the "Work-Transport" matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potential = Solution of the Poisson equation
% VFieldx_e, VFieldy_e = xy Velocity-Field components for electrons [um/ns]
% VFieldx_h, VFieldy_h = xy Velocity-Field components for holes [um/ns]
% x, y    = Axes
% Step    = Unit step of the lattice on which the field is computed [um]
% Bulk    = Bulk thickness [um]
% Radius  = Unit step of the movements and field interpolation [um]
% TauBe/h = Life-time on the backplane side [ns]
% TauSe/h = Life-time on the strip side [ns]
% ItFigIn = Figure iterator input

function [WorkTransportTotal, x, y, ItFigOut] =...
    WorkTransport(Potential,VFieldx_e,VFieldy_e,VFieldx_h,VFieldy_h,...
    x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
[mshx,mshy] = meshgrid(x,y);
queryPoints = [mshx(:),mshy(:)]';

[Ewx, Ewy] = evaluateGradient(Potential,queryPoints);

Ewx = reshape(Ewx,size(mshx));
Ewy = reshape(Ewy,size(mshy));

interp = interpolateSolution(Potential,queryPoints);
interp = reshape(interp,size(mshx));

fprintf('@@@ I''m calculating the work for electrons @@@\n');
[WTransport_e] =...
    WorkTransportMatrix(VFieldx_e,VFieldy_e,Ewx,Ewy,x,y,...
    TauBe,TauSe,Step,Bulk,Radius,-1);

fprintf('@@@ I''m calculating the work for holes @@@\n');
[WTransport_h] =...
    WorkTransportMatrix(VFieldx_h,VFieldy_h,Ewx,Ewy,x,y,...
    TauBh,TauSh,Step,Bulk,Radius,+1);

WorkTransportTotal = WTransport_e + WTransport_h;


%%%%%%%%%
% Plots %
%%%%%%%%%
[xx,yy] = meshgrid(x,y);

figure(ItFigIn);
colormap jet;
surf(xx,yy,interp,'EdgeColor','none','FaceColor','interp');
title('Weighting potential');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('Potential [V]');

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
colormap jet;
surf(xx,yy,WorkTransportTotal,'EdgeColor','none','FaceColor','interp');
title('Total Work-Transport');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('Work / q [#charges * V]');

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
colormap jet;
subplot(1,2,1);
surf(xx,yy,WTransport_e,'EdgeColor','none','FaceColor','interp');
title('Electron Work-Transport');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('Work / q [#charges * V]');
subplot(1,2,2);
surf(xx,yy,WTransport_h,'EdgeColor','none','FaceColor','interp');
title('Hole Work-Transport');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('Work / q [#charges * V]');

ItFigOut = ItFigIn + 1;
fprintf('CPU time --> %.2f [min]\n',(cputime-TStart)/60);
end
