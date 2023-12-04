%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots related to solution fo 2D Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdem           = PDE solver
% Potential      = Solution of Poisson equation
% DecomposedGeom = Decomposed geometry
% BulkStart      = Bulk start coordinate [um]
% BulkStart      = Bulk stop coordinate [um]
% VolumeHeight   = Volume height [units of bulk thickness]
% PitchX         = Pitch along X [um]
% PitchY         = Pitch along Y [um]
% BiasV          = Sensor central strip voltage [V]
% epsR           = Relative permittivity
% StripNot3D == true ? Strip : 3D
% XQ             = Coordinate for potential query along y [um]
% ItFigIn        = Figure iterator input

function [Sq, yq, ItFigOut] = Planar_Plots(pdem,Potential,DecomposedGeom,...
    BulkStart,BulkStop,VolumeHeight,PitchX,PitchY,BiasV,epsR,StripNot3D,XQ,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps0           = 8.85e-18; % Vacuum permittivity [F/um]
ReSampleFine   = 1;        % Used in order to make nice plots [um]
ReSampleCoarse = 10;       % Used in order to make nice plots [um]
ContLevel      = 40;       % Contour plot levels
MagnVector     = 0.5;      % Vector field magnification
NStrips        = 13;       % Total number of strips
NPixelsX       = 5;        % Number of pixels along X
NPixelsY       = 5;        % Number of pixels along Y


%%%%%%%%%%%%%%%%%%%%%%%
% Specific parameters %
%%%%%%%%%%%%%%%%%%%%%%%
if StripNot3D == true
    % Planar Strip
    myXlim  = [-PitchX * NStrips/2,+PitchX * NStrips/2];
    myYlim  = [0,(BulkStop - BulkStart) * VolumeHeight];
    L       = PitchY;
    xfine   = -PitchX:ReSampleFine:PitchX;
    yfine   = 0:ReSampleFine:(BulkStop - BulkStart) * VolumeHeight;
    xcoarse = -PitchX:ReSampleCoarse:PitchX;
    ycoarse = 0:ReSampleCoarse:(BulkStop - BulkStart) * VolumeHeight;
else
    % 3D Pixel
    myXlim   = [-(PitchX*NPixelsX/2+PitchX/2),PitchX*NPixelsX/2+PitchX/2];
    myYlim   = [-(PitchY*NPixelsY/2+PitchY/2),PitchY*NPixelsY/2+PitchY/2];
    L        = BulkStop - BulkStart;
    xfine    = -PitchX/2:ReSampleFine:PitchX/2;
    yfine    = -PitchY/2:ReSampleFine:PitchY/2;
    xcoarse  = -PitchX/2:ReSampleCoarse:PitchX/2;
    ycoarse  = -PitchY/2:ReSampleCoarse:PitchY/2;
end


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
interp = interpolateElectricPotential(Potential,FineQuery);
interp = reshape(interp,size(FineMeshX));

FineGrad   = interpolateElectricField(Potential,FineQuery);
CoarseGrad = interpolateElectricField(Potential,CoarseQuery);

EfieldNorm = reshape(sqrt(FineGrad.Ex.^2 + FineGrad.Ey.^2),size(FineMeshX));
EfieldNorm(isinf(EfieldNorm) | isnan(EfieldNorm)) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate capacitance %
%%%%%%%%%%%%%%%%%%%%%%%%
U = trapz(yfine,trapz(xfine,1/2 * EfieldNorm .* EfieldNorm,2));
C = eps0*epsR * 2*U / (BiasV * BiasV) / 1e-12; % Capacitance [pF/um]
fprintf('Channel capacitance --> %.4f [pF/um] --> %.2f [pF]\n',C,C*L);


%%%%%%%%%
% Plots %
%%%%%%%%%
surf(FineMeshX,FineMeshY,EfieldNorm,'FaceAlpha',0.9,'EdgeColor','none','FaceColor','interp');
hold on;
contour(FineMeshX,FineMeshY,interp,ContLevel);
quiver(CoarseMeshX(:),CoarseMeshY(:),CoarseGrad.Ex,CoarseGrad.Ey,MagnVector);
hold off;
title('Potential, gradient, and gradient magnitude');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Electric field abs. value [V/\mum]');
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate potential along a line %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if StripNot3D == true
    ItFigIn = ItFigIn + 1;
    figure(ItFigIn);
    yq = BulkStart:ReSampleFine:BulkStop;
    xq = XQ * ones(1,length(yq));
    Sq = interpolateElectricPotential(Potential,xq,yq);
    plot(yq,Sq);
    title(sprintf('Potential along z at x = %.2f um',XQ));
    xlabel('Z [\mum]');
    ylabel('Potential');
    grid on;
else
    Sq = 0;
    yq = 0;
end


ItFigOut = ItFigIn + 1;
fprintf('CPU time --> %.2f [min]\n\n',(cputime-TStart)/60);
end
