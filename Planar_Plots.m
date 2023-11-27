%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots related to solution fo 2D Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdem      = PDE solver
% Potential = Splution Poisson equation
% DecomposedGeom = Decomposed geometry
% Bulk      = Bulk thickness [um]
% BulkStart = Bulk start coordinate [um]
% BulkStart = Bulk stop coordinate [um]
% PitchX    = Pitch along X [um]
% PitchY    = Pitch along Y [um]
% BiasW     = Sensor central strip voltage [V] [1 Weighting; 0 All]
% epsR      = Relative permittivity
% XQ        = Coordinate for potential query along y [um]
% ItFigIn   = Figure iterator input

function [Sq, yq, ItFigOut] = Planar_Plots(pdem,Potential,DecomposedGeom,...
    Bulk,BulkStart,BulkStop,PitchX,PitchY,BiasW,epsR,XQ,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps0           = 8.85e-18; % Vacuum permittivity [F/um]
ReSampleFine   = 1;        % Used in order to make nice plots [um]
ReSampleCoarse = 10;       % Used in order to make nice plots [um]
ContLevel      = 40;       % Contour plot levels
MagnVector     = 1.5;      % Vector field magnification
VolumeHeight   = 2;        % Volume height [units of bulk thickness]
NStrips        = 13;       % Total number of strips


%%%%%%%%%%%%%%%%%%%%%%%
% Specific parameters %
%%%%%%%%%%%%%%%%%%%%%%%

% Planar Pixel
%myXlim  = -(PitchX*NPixelsX/2+PitchX/2),PitchX*NPixelsX/2+PitchX/2];
%myYlim  = -(PitchY*NPixelsY/2+PitchY/2),PitchY*NPixelsY/2+PitchY/2];
%L       = Bulk;
%xfine   = -PitchX:ReSampleFine:PitchX;
%yfine   = -PitchY:ReSampleFine:PitchY;
%xcoarse = -PitchX:ReSampleCoarse:PitchX;
%ycoarse = -PitchY:ReSampleCoarse:PitchY;

% Planar Strip
myXlim  = [-PitchX * NStrips/2,+PitchX * NStrips/2];
myYlim  = [0,Bulk * VolumeHeight];
L       = PitchY;
xfine   = -PitchX:ReSampleFine:PitchX;
yfine   = BulkStart:ReSampleFine:BulkStop * 3/2;
xcoarse = -PitchX:ReSampleCoarse:PitchX;
ycoarse = BulkStart:ReSampleCoarse:BulkStop * 3/2;


%%%%%%%%%
% Plots %
%%%%%%%%%
figure(ItFigIn);
subplot(1,2,1);
pdegplot(DecomposedGeom,'EdgeLabels','on','SubdomainLabels','on');
xlim(myXlim);
ylim(myYlim);
title('Geometry');
xlabel('X [\mum]');
ylabel('Z [\mum]');
subplot(1,2,2);
pdegplot(pdem);
hold on;
pdemesh(pdem);
xlim(myXlim);
ylim(myYlim);
hold off;
title('Delaunay mesh');
xlabel('X [\mum]');
ylabel('Z [\mum]');

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
subplot(1,2,1);
colormap jet;
pdeplot(pdem,'xydata',Potential.NodalSolution);
xlim(myXlim);
ylim(myYlim);
title('Potential');
xlabel('X [\mum]');
ylabel('Z [\mum]');

subplot(1,2,2);
colormap jet;


%%%%%%%%%%%%%%%%%
% Redefine mesh %
%%%%%%%%%%%%%%%%%
[FineMeshX,FineMeshY] = meshgrid(xfine,yfine);
FineQuery = [FineMeshX(:),FineMeshY(:)]';

[CoarseMeshX,CoarseMeshY] = meshgrid(xcoarse,ycoarse);
CoarseQuery = [CoarseMeshX(:),CoarseMeshY(:)]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate gradient field %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
interp = interpolateSolution(Potential,FineQuery);
interp = reshape(interp,size(FineMeshX));
[FineGradx,FineGrady]     = evaluateGradient(Potential,FineQuery);
[CoarseGradx,CoarseGrady] = evaluateGradient(Potential,CoarseQuery);
EfieldNorm = reshape(sqrt(FineGradx.^2 + FineGrady.^2),size(FineMeshX));
EfieldNorm(isinf(EfieldNorm) | isnan(EfieldNorm)) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate capacitance %
%%%%%%%%%%%%%%%%%%%%%%%%
U = trapz(yfine,trapz(xfine,1/2 * EfieldNorm .* EfieldNorm,2));
C = eps0*epsR * 2*U / (BiasW * BiasW) / 1e-12; % Capacitance [pF/um]
fprintf('Strip capacitance --> %.4f [pF/um] --> %.2f [pF]\n',C,C*L);


%%%%%%%%%
% Plots %
%%%%%%%%%
surf(FineMeshX,FineMeshY,EfieldNorm,'FaceAlpha',0.9,'EdgeColor','none','FaceColor','interp');
hold on;
contour(FineMeshX,FineMeshY,interp,ContLevel);
quiver(CoarseMeshX(:),CoarseMeshY(:),CoarseGradx,CoarseGrady,MagnVector);
hold off;
title('Potential, gradient, and gradient magnitude');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Electric field abs. value [V/\mum]');
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate potential along a line %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ItFigIn = ItFigIn + 1;
figure(ItFigIn);
yq = BulkStart:ReSampleFine:BulkStop;
xq = XQ * ones(1,length(yq));
Sq = interpolateSolution(Potential,xq,yq);
plot(yq,Sq);
title(sprintf('Potential along z at x = %.2f um',XQ));
xlabel('Z [\mum]');
ylabel('Potential');
grid on;

ItFigOut = ItFigIn + 1;
fprintf('CPU time --> %.2f [min]\n\n',(cputime-TStart)/60);
end
