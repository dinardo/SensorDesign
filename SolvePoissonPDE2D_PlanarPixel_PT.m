%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation to compute the potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk     = Bulk thickness [um]
% Pitch    = Strip pitch [um]
% BiasB    = Sensor backplane voltage [V] [0 Weighting; -V All]
% BiasS    = Sensor strip voltage [V]
% BiasG    = Bias grid voltage [V]
% epsR     = Relative permittivity of the bulk
% epsRSiO2 = Relative permittivity of the bias grid insulation (SiO2)
% rho      = Charge density in the bulk [(Coulomb/um^3) / eps0 [F/um]]
% XQ       = Coordinate for potential query along y [um]
% ItFigIn  = Figure iterator input

function [Potential, Sq, yq, ItFigOut] = SolvePoissonPDE2D_PlanarPixel_PT(Bulk,Pitch,...
    BiasB,BiasS,BiasG,epsR,epsRSiO2,rho,XQ,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReSampleFine   = 0.8; % Used in order to make nice plots [um]
ReSampleCoarse = 10;  % Used in order to make nice plots [um]
ContLevel      = 40;  % Contour plot levels
MagnVector     = 1;   % Vector field magnification
MeshMax        = 15;  % Maximum mesh edge length [um]

SHeight         = 2;  % Sensor height [units of bulk thickness]
MetalThick      = 2;  % Metalization thickness [um]
MetalWidth      = 15; % Bias grid width [um]
NStrips         = 12; % Total number of strips
InsulationThick = 2;  % SiO2 thickness [um]

%%%%%%%%%%%%%%%%%%%%
% Create PDE model %
%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m solving Poisson equation in 2D to calculate the potential @@@\n');
pdem = createpde(1);


%%%%%%%%%%%%%%%%%%%%%%
% Create 2D geometry %
%%%%%%%%%%%%%%%%%%%%%%
% Central strip
R1 = [ 3 4 -MetalWidth/2 MetalWidth/2 MetalWidth/2 -MetalWidth/2 ...
    Bulk+MetalThick+InsulationThick Bulk+MetalThick+InsulationThick Bulk+InsulationThick Bulk+InsulationThick ]';
% Positive strips
R2 = [ 3 4 1*Pitch/2 NStrips/2*Pitch NStrips/2*Pitch 1*Pitch/2 ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
% Negative strips
R3 = [ 3 4 -(1*Pitch/2) -(NStrips/2*Pitch) -(NStrips/2*Pitch) -(1*Pitch/2) ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
% Insulation
R4 = [ 3 4 -Pitch/2 Pitch/2 Pitch/2 -Pitch/2 ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
% Sensor volume
R5 = [ 3 4 -(NStrips/2*Pitch) (NStrips/2*Pitch) (NStrips/2*Pitch) -(NStrips/2*Pitch) ...
    Bulk Bulk 0 0 ]';
% Whole volume
R6 = [ 3 4 -(NStrips/2*Pitch) (NStrips/2*Pitch) (NStrips/2*Pitch) -(NStrips/2*Pitch) ...
    Bulk*SHeight Bulk*SHeight 0 0 ]';

gd = [R1,R2,R3,R4,R5,R6];
sf = 'R6-(R1+R2+R3)+R4+R5';
ns = char('R1','R2','R3','R4','R5','R6');
ns = ns';

dl = decsg(gd,sf,ns);
geometryFromEdges(pdem,dl);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply boundary conditions (only on conductors) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Central strip
applyBoundaryCondition(pdem,'edge',1,'h',1,'r',BiasG);
applyBoundaryCondition(pdem,'edge',2,'h',1,'r',BiasG);
applyBoundaryCondition(pdem,'edge',3,'h',1,'r',BiasG);
applyBoundaryCondition(pdem,'edge',7,'h',1,'r',BiasG);
% Positive strips
applyBoundaryCondition(pdem,'edge',9,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',10,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',14,'h',1,'r',0);
% Negative strips
applyBoundaryCondition(pdem,'edge',5,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',12,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',11,'h',1,'r',0);
% Top edge
applyBoundaryCondition(pdem,'edge',4,'h',1,'r',BiasS);
% Right edge sensor
applyBoundaryCondition(pdem,'edge',15,'q',0,'g',0);
% Right edge air
applyBoundaryCondition(pdem,'edge',16,'q',0,'g',0);
% Bottom edge
applyBoundaryCondition(pdem,'edge',17,'h',1,'r',BiasB);
% Left edge sensor
applyBoundaryCondition(pdem,'edge',19,'q',0,'g',0);
% Left edge air
applyBoundaryCondition(pdem,'edge',18,'q',0,'g',0);


%%%%%%%%%%%%%%%%%
% Generate mesh %
%%%%%%%%%%%%%%%%%
msh = generateMesh(pdem,'Hmax',MeshMax,'Jiggle','mean',...
    'GeometricOrder','quadratic','MesherVersion','R2013a');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%
specifyCoefficients(pdem,'m',0,'d',0,'c',1,       'a',0,'f',0,  'face',1); % Air
specifyCoefficients(pdem,'m',0,'d',0,'c',epsR,    'a',0,'f',rho,'face',2); % Sensor
specifyCoefficients(pdem,'m',0,'d',0,'c',epsRSiO2,'a',0,'f',0,  'face',3); % Insulator
Potential = solvepde(pdem);


%%%%%%%%%
% Plots %
%%%%%%%%%
figure(ItFigIn);
subplot(1,2,1);
pdegplot(dl,'EdgeLabels','on','SubdomainLabels','on');
xlim([-Pitch * NStrips/2,+Pitch * NStrips/2]);
ylim([0,Bulk * SHeight]);
title('Geometry');
xlabel('X [\mum]');
ylabel('Z [\mum]');
subplot(1,2,2);
pdegplot(pdem);
hold on;
pdemesh(pdem);
xlim([-Pitch * NStrips/2,+Pitch * NStrips/2]);
ylim([0,Bulk * SHeight]);
hold off;
title('Delaunay mesh');
xlabel('X [\mum]');
ylabel('Z [\mum]');

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
subplot(1,2,1);
colormap jet;
pdeplot(pdem,'xydata',Potential.NodalSolution);
xlim([-Pitch * NStrips/2,+Pitch * NStrips/2]);
ylim([0,Bulk * SHeight]);
title('Potential');
xlabel('X [\mum]');
ylabel('Z [\mum]');


%%%%%%%%%%%%%%%%%
% Redefine mesh %
%%%%%%%%%%%%%%%%%
xfine = -Pitch:ReSampleFine:Pitch;
yfine = 0:ReSampleFine:Bulk * 3/2;
[FineMeshX,FineMeshY] = meshgrid(xfine,yfine);
FineQuery = [FineMeshX(:),FineMeshY(:)]';

xcoarse = -Pitch:ReSampleCoarse:Pitch;
ycoarse = 0:ReSampleCoarse:Bulk * 3/2;
[CoarseMeshX,CoarseMeshY] = meshgrid(xcoarse,ycoarse);
CoarseQuery = [CoarseMeshX(:),CoarseMeshY(:)]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recompute solution on a different mesh %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interp = interpolateSolution(Potential,FineQuery);
interp = reshape(interp,size(FineMeshX));


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate gradient field %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[gradx,grady] = evaluateGradient(Potential,FineQuery);
contour(FineMeshX,FineMeshY,interp,ContLevel);

subplot(1,2,2);
colormap jet;
hold on;
quiver(FineMeshX(:),FineMeshY(:),gradx,grady,MagnVector);
hold off;
title('Potential and its gradient');
xlabel('X [\mum]');
ylabel('Z [\mum]');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate gradient magnitude %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[gradx,grady] = evaluateGradient(Potential,FineQuery);
gradx = reshape(gradx,length(xfine),length(yfine));
grady = reshape(grady,length(xfine),length(yfine));
EfieldNorm = sqrt(gradx.^2 + grady.^2);
xx = reshape(FineQuery(1,:),length(xfine),length(yfine));
yy = reshape(FineQuery(2,:),length(xfine),length(yfine));

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
surf(xx,yy,EfieldNorm,'EdgeColor','none','FaceColor','interp');
title('Field magnitude');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('|E| [V/\mum]');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate potential along a line %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ItFigIn = ItFigIn + 1;
figure(ItFigIn);
yq = 0:ReSampleFine:Bulk;
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
