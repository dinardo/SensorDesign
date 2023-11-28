%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation to compute the potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk    = Bulk thickness [um]
% PitchX  = Pitch along X [um]
% BiasB   = Sensor backplane voltage [V] [0 Weighting; -V All]
% BiasW   = Sensor central strip voltage [V] [1 Weighting; 0 All]
% epsR    = Relative permittivity
% rho     = Charge density in the bulk [(Coulomb/um^3) / eps0 [F/um]]

function [pdem,Potential,DecomposedGeom,BulkStart,BulkStop] = StripPlanarDouble_SolvePoisson2D(...
    Bulk,PitchX,BiasB,BiasW,epsR,rho)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
MeshMax      = 5;           % Maximum mesh edge length [um]
VolumeHeight = 3;           % Volume height [units of bulk thickness]
MetalThick   = 5;           % Metalization thickness [um]
MetalWidth   = PitchX-20;   % Metalization width [um]
BulkStart    = Bulk/2;      % Bulk start coordinate [um]
BulkStop     = Bulk/2+Bulk; % Bulk stop coordinate [um]


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
R2 = [ 3 4 -MetalWidth/2 MetalWidth/2 MetalWidth/2 -MetalWidth/2 ...
    BulkStart-MetalThick BulkStart-MetalThick BulkStart BulkStart ]';
% Positive up-strips
R3 = [ 3 4 1*PitchX-MetalWidth/2 1*PitchX+MetalWidth/2 1*PitchX+MetalWidth/2 1*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R4 = [ 3 4 2*PitchX-MetalWidth/2 2*PitchX+MetalWidth/2 2*PitchX+MetalWidth/2 2*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R5 = [ 3 4 3*PitchX-MetalWidth/2 3*PitchX+MetalWidth/2 3*PitchX+MetalWidth/2 3*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R6 = [ 3 4 4*PitchX-MetalWidth/2 4*PitchX+MetalWidth/2 4*PitchX+MetalWidth/2 4*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R7 = [ 3 4 5*PitchX-MetalWidth/2 5*PitchX+MetalWidth/2 5*PitchX+MetalWidth/2 5*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R8 = [ 3 4 6*PitchX-MetalWidth/2 6*PitchX+MetalWidth/2 6*PitchX+MetalWidth/2 6*PitchX-MetalWidth/2 ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
% Negative up-strips
R9 = [ 3 4 -(1*PitchX-MetalWidth/2) -(1*PitchX+MetalWidth/2) -(1*PitchX+MetalWidth/2) -(1*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R10 = [ 3 4 -(2*PitchX-MetalWidth/2) -(2*PitchX+MetalWidth/2) -(2*PitchX+MetalWidth/2) -(2*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R11 = [ 3 4 -(3*PitchX-MetalWidth/2) -(3*PitchX+MetalWidth/2) -(3*PitchX+MetalWidth/2) -(3*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R12 = [ 3 4 -(4*PitchX-MetalWidth/2) -(4*PitchX+MetalWidth/2) -(4*PitchX+MetalWidth/2) -(4*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R13 = [ 3 4 -(5*PitchX-MetalWidth/2) -(5*PitchX+MetalWidth/2) -(5*PitchX+MetalWidth/2) -(5*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
R14 = [ 3 4 -(6*PitchX-MetalWidth/2) -(6*PitchX+MetalWidth/2) -(6*PitchX+MetalWidth/2) -(6*PitchX-MetalWidth/2) ...
    BulkStop+MetalThick BulkStop+MetalThick BulkStop BulkStop ]';
% Positive down-strips
R15 = [ 3 4 1*PitchX-MetalWidth/2 1*PitchX+MetalWidth/2 1*PitchX+MetalWidth/2 1*PitchX-MetalWidth/2 ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
R16 = [ 3 4 2*PitchX-MetalWidth/2 2*PitchX+MetalWidth/2 2*PitchX+MetalWidth/2 2*PitchX-MetalWidth/2 ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
R17 = [ 3 4 3*PitchX-MetalWidth/2 3*PitchX+MetalWidth/2 3*PitchX+MetalWidth/2 3*PitchX-MetalWidth/2 ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
R18 = [ 3 4 4*PitchX-MetalWidth/2 4*PitchX+MetalWidth/2 4*PitchX+MetalWidth/2 4*PitchX-MetalWidth/2 ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
R19 = [ 3 4 5*PitchX-MetalWidth/2 5*PitchX+MetalWidth/2 5*PitchX+MetalWidth/2 5*PitchX-MetalWidth/2 ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
R20 = [ 3 4 6*PitchX-MetalWidth/2 6*PitchX+MetalWidth/2 6*PitchX+MetalWidth/2 6*PitchX-MetalWidth/2 ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
% Negative down-strips
R21 = [ 3 4 -(1*PitchX-MetalWidth/2) -(1*PitchX+MetalWidth/2) -(1*PitchX+MetalWidth/2) -(1*PitchX-MetalWidth/2) ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
R22 = [ 3 4 -(2*PitchX-MetalWidth/2) -(2*PitchX+MetalWidth/2) -(2*PitchX+MetalWidth/2) -(2*PitchX-MetalWidth/2) ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
R23 = [ 3 4 -(3*PitchX-MetalWidth/2) -(3*PitchX+MetalWidth/2) -(3*PitchX+MetalWidth/2) -(3*PitchX-MetalWidth/2) ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
R24 = [ 3 4 -(4*PitchX-MetalWidth/2) -(4*PitchX+MetalWidth/2) -(4*PitchX+MetalWidth/2) -(4*PitchX-MetalWidth/2) ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
R25 = [ 3 4 -(5*PitchX-MetalWidth/2) -(5*PitchX+MetalWidth/2) -(5*PitchX+MetalWidth/2) -(5*PitchX-MetalWidth/2) ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
R26 = [ 3 4 -(6*PitchX-MetalWidth/2) -(6*PitchX+MetalWidth/2) -(6*PitchX+MetalWidth/2) -(6*PitchX-MetalWidth/2) ...
    BulkStart BulkStart BulkStart-MetalThick BulkStart-MetalThick ]';
% Sensor volume
R27 = [ 3 4 -(6*PitchX+PitchX/2) (6*PitchX+PitchX/2) (6*PitchX+PitchX/2) -(6*PitchX+PitchX/2) ...
    BulkStop BulkStop BulkStart BulkStart ]';
% Whole volume
R28 = [ 3 4 -(6*PitchX+PitchX/2) (6*PitchX+PitchX/2) (6*PitchX+PitchX/2) -(6*PitchX+PitchX/2) ...
    Bulk*VolumeHeight Bulk*VolumeHeight 0 0 ]';

gd = [R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,R16,R17,R18,R19,R20,R21,R22,R23,R24,R25,R26,R27,R28];
sf = 'R28 - (R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13+R14+R15+R16+R17+R18+R19+R20+R21+R22+R23+R24+R25+R26) + R27';
ns = char('R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14','R15','R16','R17','R18','R19','R20','R21','R22','R23','R24','R25','R26','R27','R28');
ns = ns';

DecomposedGeom = decsg(gd,sf,ns);
geometryFromEdges(pdem,DecomposedGeom);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply boundary conditions (only on conductors) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Central up-strip
applyBoundaryCondition(pdem,'dirichlet','edge',4,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'dirichlet','edge',6,'h',1,'r',BiasW);
% Positive up-strips
applyBoundaryCondition(pdem,'dirichlet','edge',23,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',21,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',27,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',25,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',31,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',29,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',35,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',33,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',39,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',37,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',43,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',41,'h',1,'r',0);
% Negative up-strips
applyBoundaryCondition(pdem,'dirichlet','edge',47,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',45,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',51,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',49,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',55,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',53,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',59,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',57,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',63,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',61,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',80,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',65,'h',1,'r',0);
% Top all up-strips
applyBoundaryCondition(pdem,'dirichlet','edge',14,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',15,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',16,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',17,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',18,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',19,'h',1,'r',0);

applyBoundaryCondition(pdem,'dirichlet','edge',13,'h',1,'r',BiasW);

applyBoundaryCondition(pdem,'dirichlet','edge',12,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',11,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',10,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',9,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',8,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',7,'h',1,'r',0);
% Bottom all up-strips
applyBoundaryCondition(pdem,'dirichlet','edge',96,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',98,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',100,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',102,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',104,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',106,'h',1,'r',0);

applyBoundaryCondition(pdem,'dirichlet','edge',94,'h',1,'r',BiasW);

applyBoundaryCondition(pdem,'dirichlet','edge',92,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',90,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',88,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',86,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',84,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',82,'h',1,'r',0);
% Top edge
applyBoundaryCondition(pdem,'dirichlet','edge',1,'h',1,'r',0);
% Right edge air
applyBoundaryCondition(pdem,'neumann','edge',137,'q',0,'g',0);
% Right edge sensor
applyBoundaryCondition(pdem,'neumann','edge',136,'q',0,'g',0);
% Right edge air
applyBoundaryCondition(pdem,'neumann','edge',135,'q',0,'g',0);
% Bottom edge
applyBoundaryCondition(pdem,'dirichlet','edge',2,'h',1,'r',0);
% Left edge air
applyBoundaryCondition(pdem,'neumann','edge',140,'q',0,'g',0);
% Left edge sensor
applyBoundaryCondition(pdem,'neumann','edge',139,'q',0,'g',0);
% Left edge air
applyBoundaryCondition(pdem,'neumann','edge',138,'q',0,'g',0);
% Central down-strip
applyBoundaryCondition(pdem,'dirichlet','edge',3,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',5,'h',1,'r',BiasB);
% Positive down-strips
applyBoundaryCondition(pdem,'dirichlet','edge',22,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',20,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',26,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',24,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',30,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',28,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',34,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',32,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',38,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',36,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',42,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',40,'h',1,'r',BiasB);
% Negative down-strips
applyBoundaryCondition(pdem,'dirichlet','edge',46,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',44,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',50,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',48,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',54,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',52,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',58,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',56,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',62,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',60,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',79,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',64,'h',1,'r',BiasB);
% Top all down-strips
applyBoundaryCondition(pdem,'dirichlet','edge',123,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',125,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',127,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',129,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',131,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',133,'h',1,'r',BiasB);

applyBoundaryCondition(pdem,'dirichlet','edge',121,'h',1,'r',BiasB);

applyBoundaryCondition(pdem,'dirichlet','edge',119,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',117,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',115,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',113,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',111,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',109,'h',1,'r',BiasB);
% Bottom all down-strips
applyBoundaryCondition(pdem,'dirichlet','edge',73,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',74,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',75,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',76,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',77,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',78,'h',1,'r',BiasB);

applyBoundaryCondition(pdem,'dirichlet','edge',72,'h',1,'r',BiasB);

applyBoundaryCondition(pdem,'dirichlet','edge',71,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',70,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',69,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',68,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',67,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'dirichlet','edge',66,'h',1,'r',BiasB);


%%%%%%%%%%%%%%%%%
% Generate mesh %
%%%%%%%%%%%%%%%%%
msh = generateMesh(pdem,'Hmax',MeshMax,'GeometricOrder','quadratic');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%
specifyCoefficients(pdem,'m',0,'d',0,'c',1,   'a',0,'f',0,  'face',1); % Air
specifyCoefficients(pdem,'m',0,'d',0,'c',epsR,'a',0,'f',rho,'face',3); % Sensor
specifyCoefficients(pdem,'m',0,'d',0,'c',1,   'a',0,'f',0,  'face',2); % Air
Potential = solvepde(pdem);


fprintf('CPU time --> %.2f [min]\n\n',(cputime-TStart)/60);
end
