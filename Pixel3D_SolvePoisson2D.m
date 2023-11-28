%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation to compute the potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PitchX = Pitch along X [um]
% PitchY = Pitch along Y [um]
% BiasB  = Sensor backplane voltage [V] [0 Weighting; -V All]
% BiasW  = Sensor central pixel voltage [V] [1 Weighting; 0 All]
% epsR   = Relative permittivity
% rho    = Charge density in the bulk [(Coulomb/um^3) / eps0 [F/um]]

function [pdem,Potential,DecomposedGeom] = Pixel3D_SolvePoisson2D(...
    PitchX,PitchY,BiasB,BiasW,epsR,rho)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
MeshMax  = 5;   % Maximum mesh edge length [um]
Radius   = 2.5; % Column radius [um]
NPixelsX = 5;   % Number of pixels along X
NPixelsY = 5;   % Number of pixels along Y


%%%%%%%%%%%%%%%%%%%%
% Create PDE model %
%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m solving Poisson equation in 2D to calculate the potential @@@\n');
pdem = createpde(1);


%%%%%%%%%%%%%%%%%%%%%%
% Create 2D geometry %
%%%%%%%%%%%%%%%%%%%%%%
% Whole volume
R1 = [ 3 4 -(PitchX*NPixelsX/2+PitchX/2) PitchX*NPixelsX/2+PitchX/2 PitchX*NPixelsX/2+PitchX/2 -(PitchX*NPixelsX/2+PitchX/2)...
    PitchY*NPixelsY/2+PitchY/2 PitchY*NPixelsY/2+PitchY/2 -(PitchY*NPixelsY/2+PitchY/2) -(PitchY*NPixelsY/2+PitchY/2) ]';

%%%%%%%%%%%%%%%%%%%
% Junction colums %
%%%%%%%%%%%%%%%%%%%
Cs11 = [ 1  PitchX*0  PitchY*0 Radius 0 0 0 0 0 0 ]';
Cs12 = [ 1  PitchX*2 -PitchY*1 Radius 0 0 0 0 0 0 ]';
Cs13 = [ 1  PitchX*2  PitchY*0 Radius 0 0 0 0 0 0 ]';
Cs14 = [ 1  PitchX*2  PitchY*1 Radius 0 0 0 0 0 0 ]';
Cs15 = [ 1  PitchX*2  PitchY*2 Radius 0 0 0 0 0 0 ]';

Cs21 = [ 1  PitchX*1 -PitchY*2 Radius 0 0 0 0 0 0 ]';
Cs22 = [ 1  PitchX*1 -PitchY*1 Radius 0 0 0 0 0 0 ]';
Cs23 = [ 1  PitchX*1  PitchY*0 Radius 0 0 0 0 0 0 ]';
Cs24 = [ 1  PitchX*1  PitchY*1 Radius 0 0 0 0 0 0 ]';
Cs25 = [ 1  PitchX*1  PitchY*2 Radius 0 0 0 0 0 0 ]';

Cs31 = [ 1  PitchX*0 -PitchY*2 Radius 0 0 0 0 0 0 ]';
Cs32 = [ 1  PitchX*0 -PitchY*1 Radius 0 0 0 0 0 0 ]';
Cs33 = [ 1  PitchX*2 -PitchY*2 Radius 0 0 0 0 0 0 ]'; % Scambio con Cs11
Cs34 = [ 1  PitchX*0  PitchY*1 Radius 0 0 0 0 0 0 ]';
Cs35 = [ 1  PitchX*0  PitchY*2 Radius 0 0 0 0 0 0 ]';

Cs41 = [ 1 -PitchX*1 -PitchY*2 Radius 0 0 0 0 0 0 ]';
Cs42 = [ 1 -PitchX*1 -PitchY*1 Radius 0 0 0 0 0 0 ]';
Cs43 = [ 1 -PitchX*1  PitchY*0 Radius 0 0 0 0 0 0 ]';
Cs44 = [ 1 -PitchX*1  PitchY*1 Radius 0 0 0 0 0 0 ]';
Cs45 = [ 1 -PitchX*1  PitchY*2 Radius 0 0 0 0 0 0 ]';

Cs51 = [ 1 -PitchX*2 -PitchY*2 Radius 0 0 0 0 0 0 ]';
Cs52 = [ 1 -PitchX*2 -PitchY*1 Radius 0 0 0 0 0 0 ]';
Cs53 = [ 1 -PitchX*2  PitchY*0 Radius 0 0 0 0 0 0 ]';
Cs54 = [ 1 -PitchX*2  PitchY*1 Radius 0 0 0 0 0 0 ]';
Cs55 = [ 1 -PitchX*2  PitchY*2 Radius 0 0 0 0 0 0 ]';

%%%%%%%%%%%%%%%
% Bias colums %
%%%%%%%%%%%%%%%
Cb11 = [ 1  PitchX*2+PitchX/2 -PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb12 = [ 1  PitchX*2+PitchX/2 -PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb13 = [ 1  PitchX*2+PitchX/2  PitchY*0+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb14 = [ 1  PitchX*2+PitchX/2  PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb15 = [ 1  PitchX*2+PitchX/2  PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';

Cb21 = [ 1  PitchX*1+PitchX/2 -PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb22 = [ 1  PitchX*1+PitchX/2 -PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb23 = [ 1  PitchX*1+PitchX/2  PitchY*0+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb24 = [ 1  PitchX*1+PitchX/2  PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb25 = [ 1  PitchX*1+PitchX/2  PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';

Cb31 = [ 1  PitchX*0+PitchX/2 -PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb32 = [ 1  PitchX*0+PitchX/2 -PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb33 = [ 1  PitchX*0+PitchX/2  PitchY*0+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb34 = [ 1  PitchX*0+PitchX/2  PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb35 = [ 1  PitchX*0+PitchX/2  PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';

Cb41 = [ 1 -PitchX*1+PitchX/2 -PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb42 = [ 1 -PitchX*1+PitchX/2 -PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb43 = [ 1 -PitchX*1+PitchX/2  PitchY*0+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb44 = [ 1 -PitchX*1+PitchX/2  PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb45 = [ 1 -PitchX*1+PitchX/2  PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';

Cb51 = [ 1 -PitchX*2+PitchX/2 -PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb52 = [ 1 -PitchX*2+PitchX/2 -PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb53 = [ 1 -PitchX*2+PitchX/2  PitchY*0+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb54 = [ 1 -PitchX*2+PitchX/2  PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb55 = [ 1 -PitchX*2+PitchX/2  PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';

%%%%%%%%%%%%%%%%%%%%%
% Extra bias colums %
%%%%%%%%%%%%%%%%%%%%%
Cb61 = [ 1 -PitchX*3+PitchX/2 -PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb62 = [ 1 -PitchX*3+PitchX/2 -PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb63 = [ 1 -PitchX*3+PitchX/2  PitchY*0+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb64 = [ 1 -PitchX*3+PitchX/2  PitchY*1+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb65 = [ 1 -PitchX*3+PitchX/2  PitchY*2+PitchY/2 Radius 0 0 0 0 0 0 ]';

Cb16 = [ 1  PitchX*2+PitchX/2 -PitchY*3+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb26 = [ 1  PitchX*1+PitchX/2 -PitchY*3+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb36 = [ 1  PitchX*0+PitchX/2 -PitchY*3+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb46 = [ 1 -PitchX*1+PitchX/2 -PitchY*3+PitchY/2 Radius 0 0 0 0 0 0 ]';
Cb56 = [ 1 -PitchX*2+PitchX/2 -PitchY*3+PitchY/2 Radius 0 0 0 0 0 0 ]';

Cb66 = [ 1 -PitchX*3+PitchX/2 -PitchY*3+PitchY/2 Radius 0 0 0 0 0 0 ]';

gd = [R1,...
    Cs11,Cs12,Cs13,Cs14,Cs15,...
    Cs21,Cs22,Cs23,Cs24,Cs25,...
    Cs31,Cs32,Cs33,Cs34,Cs35,...
    Cs41,Cs42,Cs43,Cs44,Cs45,...
    Cs51,Cs52,Cs53,Cs54,Cs55,...
    Cb11,Cb12,Cb13,Cb14,Cb15,...
    Cb21,Cb22,Cb23,Cb24,Cb25,...
    Cb31,Cb32,Cb33,Cb34,Cb35,...
    Cb41,Cb42,Cb43,Cb44,Cb45,...
    Cb51,Cb52,Cb53,Cb54,Cb55,...
    Cb61,Cb62,Cb63,Cb64,Cb65,...
    Cb16,Cb26,Cb36,Cb46,Cb56,...
    Cb66];

sf = 'R1-(';
ns = char('R1  ');
for i = 1:NPixelsX
    for j = 1:NPixelsY
        sf = strcat(sf,sprintf('Cs%d%d+',i,j));
        ns = [ns; sprintf('Cs%d%d',i,j)];
    end
end
for i = 1:NPixelsX+1
    for j = 1:NPixelsY+1
        sf = strcat(sf,sprintf('Cb%d%d+',i,j));
        ns = [ns; sprintf('Cb%d%d',i,j)];
    end
end
sf = sf(1:end-1);
sf = strcat(sf,')');
ns = ns';

DecomposedGeom = decsg(gd,sf,ns);
geometryFromEdges(pdem,DecomposedGeom);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply boundary conditions (only on conductors) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary of the domain
applyBoundaryCondition(pdem,'neumann','edge',1:4,'q',0,'g',0);
% Central pixel signal columns
applyBoundaryCondition(pdem,'dirichlet','edge',5:8,'h',1,'r',BiasW);
% Other pixels signal columns
applyBoundaryCondition(pdem,'dirichlet','edge',9:12 + 4*(NPixelsX*NPixelsY-2),'h',1,'r',0);
% All pixels bias columns
applyBoundaryCondition(pdem,'dirichlet','edge',13 + 4*(NPixelsX*NPixelsY-2):pdem.Geometry.NumEdges,'h',1,'r',BiasB);


%%%%%%%%%%%%%%%%%
% Generate mesh %
%%%%%%%%%%%%%%%%%
msh = generateMesh(pdem,'Hmax',MeshMax,'GeometricOrder','quadratic');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%
specifyCoefficients(pdem,'m',0,'d',0,'c',epsR,'a',0,'f',rho,'face',1);
Potential = solvepde(pdem);


fprintf('CPU time --> %.2f [min]\n\n',(cputime-TStart)/60);
end
