%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation to compute the potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk    = Bulk thickness [um]
% PitchX  = Pitch along X [um]
% BiasB   = Sensor backplane voltage [V] [0 Weighting; -V All]
% BiasW   = Sensor central strip voltage [V] [1 Weighting; 0 All]
% epsR    = Relative permittivity
% rho     = Charge density in the bulk [(Coulomb/um^3) / eps0 [F/um]]

function [pdem,Potential,DecomposedGeom,BulkStart,BulkStop] = StripPlanarSiO2_SolvePoisson2D(...
    Bulk,PitchX,BiasB,BiasW,epsR,rho)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
MeshMax      = 5;         % Maximum mesh edge length [um]
VolumeHeight = 3;         % Volume height [units of bulk thickness]
MetalThick   = 5;         % Metalization thickness [um]
MetalWidth   = PitchX-20; % Metalization width [um]
BulkStart    = 0;         % Bulk start coordinate [um]
BulkStop     = Bulk;      % Bulk stop coordinate [um]
epsSiO2      = 3.9;       % Relative permittivity SiO2
SiO2Thick    = Bulk-40;   % SiO2 thickness  [um]


%%%%%%%%%%%%%%%%%%%%
% Create PDE model %
%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m solving Poisson equation in 2D to calculate the potential @@@\n');
pdem = createpde(1);


%%%%%%%%%%%%%%%%%%%%%%
% Create 2D geometry %
%%%%%%%%%%%%%%%%%%%%%%
% Central strips
R1 = [ 3 4 -MetalWidth/2 MetalWidth/2 MetalWidth/2 -MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
% Positive up-strips
R2 = [ 3 4 1*PitchX-MetalWidth/2 1*PitchX+MetalWidth/2 1*PitchX+MetalWidth/2 1*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R3 = [ 3 4 2*PitchX-MetalWidth/2 2*PitchX+MetalWidth/2 2*PitchX+MetalWidth/2 2*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R4 = [ 3 4 3*PitchX-MetalWidth/2 3*PitchX+MetalWidth/2 3*PitchX+MetalWidth/2 3*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R5 = [ 3 4 4*PitchX-MetalWidth/2 4*PitchX+MetalWidth/2 4*PitchX+MetalWidth/2 4*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R6 = [ 3 4 5*PitchX-MetalWidth/2 5*PitchX+MetalWidth/2 5*PitchX+MetalWidth/2 5*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R7 = [ 3 4 6*PitchX-MetalWidth/2 6*PitchX+MetalWidth/2 6*PitchX+MetalWidth/2 6*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
% Negative up-strips
R8 = [ 3 4 -(1*PitchX-MetalWidth/2) -(1*PitchX+MetalWidth/2) -(1*PitchX+MetalWidth/2) -(1*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R9 = [ 3 4 -(2*PitchX-MetalWidth/2) -(2*PitchX+MetalWidth/2) -(2*PitchX+MetalWidth/2) -(2*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R10 = [ 3 4 -(3*PitchX-MetalWidth/2) -(3*PitchX+MetalWidth/2) -(3*PitchX+MetalWidth/2) -(3*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R11 = [ 3 4 -(4*PitchX-MetalWidth/2) -(4*PitchX+MetalWidth/2) -(4*PitchX+MetalWidth/2) -(4*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R12 = [ 3 4 -(5*PitchX-MetalWidth/2) -(5*PitchX+MetalWidth/2) -(5*PitchX+MetalWidth/2) -(5*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R13 = [ 3 4 -(6*PitchX-MetalWidth/2) -(6*PitchX+MetalWidth/2) -(6*PitchX+MetalWidth/2) -(6*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
% Sensor volume
R14 = [ 3 4 -(6*PitchX+PitchX/2) (6*PitchX+PitchX/2) (6*PitchX+PitchX/2) -(6*PitchX+PitchX/2) ...
    BulkStop BulkStop BulkStart BulkStart ]';
% Whole volume
R15 = [ 3 4 -(6*PitchX+PitchX/2) (6*PitchX+PitchX/2) (6*PitchX+PitchX/2) -(6*PitchX+PitchX/2) ...
    Bulk*VolumeHeight Bulk*VolumeHeight 0 0 ]';

BetweenS = PitchX - MetalWidth;

% SiO2 between positive strips
R16 = [ 3 4 1*PitchX/2-BetweenS/2 1*PitchX/2-BetweenS/2+BetweenS 1*PitchX/2-BetweenS/2+BetweenS 1*PitchX/2-BetweenS/2 ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';
R17 = [ 3 4 3*PitchX/2-BetweenS/2 3*PitchX/2-BetweenS/2+BetweenS 3*PitchX/2-BetweenS/2+BetweenS 3*PitchX/2-BetweenS/2 ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';
R18 = [ 3 4 5*PitchX/2-BetweenS/2 5*PitchX/2-BetweenS/2+BetweenS 5*PitchX/2-BetweenS/2+BetweenS 5*PitchX/2-BetweenS/2 ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';
R19 = [ 3 4 7*PitchX/2-BetweenS/2 7*PitchX/2-BetweenS/2+BetweenS 7*PitchX/2-BetweenS/2+BetweenS 7*PitchX/2-BetweenS/2 ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';
R20 = [ 3 4 9*PitchX/2-BetweenS/2 9*PitchX/2-BetweenS/2+BetweenS 9*PitchX/2-BetweenS/2+BetweenS 9*PitchX/2-BetweenS/2 ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';
R21 = [ 3 4 11*PitchX/2-BetweenS/2 11*PitchX/2-BetweenS/2+BetweenS 11*PitchX/2-BetweenS/2+BetweenS 11*PitchX/2-BetweenS/2 ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';
% SiO2 between negative strips
R22 = [ 3 4 -(1*PitchX/2-BetweenS/2) -(1*PitchX/2-BetweenS/2+BetweenS) -(1*PitchX/2-BetweenS/2+BetweenS) -(1*PitchX/2-BetweenS/2) ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';
R23 = [ 3 4 -(3*PitchX/2-BetweenS/2) -(3*PitchX/2-BetweenS/2+BetweenS) -(3*PitchX/2-BetweenS/2+BetweenS) -(3*PitchX/2-BetweenS/2) ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';
R24 = [ 3 4 -(5*PitchX/2-BetweenS/2) -(5*PitchX/2-BetweenS/2+BetweenS) -(5*PitchX/2-BetweenS/2+BetweenS) -(5*PitchX/2-BetweenS/2) ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';
R25 = [ 3 4 -(7*PitchX/2-BetweenS/2) -(7*PitchX/2-BetweenS/2+BetweenS) -(7*PitchX/2-BetweenS/2+BetweenS) -(7*PitchX/2-BetweenS/2) ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';
R26 = [ 3 4 -(9*PitchX/2-BetweenS/2) -(9*PitchX/2-BetweenS/2+BetweenS) -(9*PitchX/2-BetweenS/2+BetweenS) -(9*PitchX/2-BetweenS/2) ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';
R27 = [ 3 4 -(11*PitchX/2-BetweenS/2) -(11*PitchX/2-BetweenS/2+BetweenS) -(11*PitchX/2-BetweenS/2+BetweenS) -(11*PitchX/2-BetweenS/2) ...
    BulkStop BulkStop BulkStop-SiO2Thick BulkStop-SiO2Thick ]';

gd = [R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,R16,R17,R18,R19,R20,R21,R22,R23,R24,R25,R26,R27];
sf = 'R15 - (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13) + (R14-(R16+R17+R18+R19+R20+R21+R22+R23+R24+R25+R26+R27))';
ns = char('R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14','R15','R16','R17','R18','R19','R20','R21','R22','R23','R24','R25','R26','R27');
ns = ns';

DecomposedGeom = decsg(gd,sf,ns);
geometryFromEdges(pdem,DecomposedGeom);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply boundary conditions (only on conductors) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Central strip
applyBoundaryCondition(pdem,'dirichlet','edge',25,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'dirichlet','edge',49,'h',1,'r',BiasW);
% Positive strips
applyBoundaryCondition(pdem,'dirichlet','edge',23,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',29,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',27,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',33,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',31,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',37,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',35,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',41,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',39,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',45,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',43,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',1,'h',1,'r',0);
% Negative strips
applyBoundaryCondition(pdem,'dirichlet','edge',47,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',53,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',51,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',57,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',55,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',61,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',59,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',65,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',63,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',108,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',94,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',2,'h',1,'r',0);
% Top all strips
applyBoundaryCondition(pdem,'dirichlet','edge',11,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',12,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',13,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',14,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',15,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',16,'h',1,'r',0);

applyBoundaryCondition(pdem,'dirichlet','edge',10,'h',1,'r',BiasW);

applyBoundaryCondition(pdem,'dirichlet','edge',9,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',8,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',7,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',6,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',5,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',4,'h',1,'r',0);
% Bottom all strips
applyBoundaryCondition(pdem,'dirichlet','edge',81,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',83,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',85,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',87,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',89,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',91,'h',1,'r',0);

applyBoundaryCondition(pdem,'dirichlet','edge',79,'h',1,'r',BiasW);

applyBoundaryCondition(pdem,'dirichlet','edge',77,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',75,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',73,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',71,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',69,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',67,'h',1,'r',0);
% Top edge
applyBoundaryCondition(pdem,'dirichlet','edge',3,'h',1,'r',0);
% Right edge sensor
applyBoundaryCondition(pdem,'neumann','edge',17,'q',0,'g',0);
% Right edge air
applyBoundaryCondition(pdem,'neumann','edge',18,'q',0,'g',0);
% Bottom edge
applyBoundaryCondition(pdem,'dirichlet','edge',19,'h',1,'r',BiasB);
% Left edge sensor
applyBoundaryCondition(pdem,'neumann','edge',20,'q',0,'g',0);
% Left edge air
applyBoundaryCondition(pdem,'neumann','edge',21,'q',0,'g',0);


%%%%%%%%%%%%%%%%%
% Generate mesh %
%%%%%%%%%%%%%%%%%
msh = generateMesh(pdem,'Hmax',MeshMax,'GeometricOrder','quadratic');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%
specifyCoefficients(pdem,'m',0,'d',0,'c',1,      'a',0,'f',0,  'face',1); % Air
specifyCoefficients(pdem,'m',0,'d',0,'c',epsR,   'a',0,'f',rho,'face',2); % Sensor
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',3); %SiO2
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',4); %SiO2
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',5); %SiO2
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',6); %SiO2
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',7); %SiO2
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',8); %SiO2
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',9); %SiO2
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',10); %SiO2
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',11); %SiO2
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',12); %SiO2
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',13); %SiO2
specifyCoefficients(pdem,'m',0,'d',0,'c',epsSiO2,'a',0,'f',0,  'face',14); %SiO2
Potential = solvepde(pdem);


fprintf('CPU time --> %.2f [min]\n\n',(cputime-TStart)/60);
end
