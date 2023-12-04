%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation to compute the potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk   = Bulk thickness [um]
% PitchX = Pitch along X [um]
% BiasV  = Sensor backplane voltage [V] == 0 ? compute weighting field
% epsR   = Relative permittivity
% rho    = Charge density in the bulk [(Coulomb/um^3)]

function [pdem,Potential,DecomposedGeom,BulkStart,BulkStop,VolumeHeight] =...
    StripPlanarDouble_SolvePoisson2D(Bulk,PitchX,BiasV,epsR,rho)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps0         = 8.85e-18;    % Vacuum permittivity [F/um]
VolumeHeight = 3;           % Volume height [units of bulk thickness]
MeshMax      = 5;           % Maximum mesh edge length [um]
MetalThick   = 5;           % Metalization thickness [um]
MetalWidth   = PitchX-20;   % Metalization width [um]
BulkStart    = Bulk/2;      % Bulk start coordinate [um]
BulkStop     = Bulk/2+Bulk; % Bulk stop coordinate [um]
BiasW        = 0;           % Bias to compute weighting potential
if BiasV == 0
    BiasW = 1;
end


%%%%%%%%%%%%%%%%%%%%
% Create PDE model %
%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m solving Poisson equation in 2D to calculate the potential @@@\n');
pdem = createpde('electromagnetic','electrostatic');


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
electromagneticBC(pdem,'Voltage',BiasW,'edge',4);
electromagneticBC(pdem,'Voltage',BiasW,'edge',6);
% Positive up-strips
electromagneticBC(pdem,'Voltage',0,'edge',23);
electromagneticBC(pdem,'Voltage',0,'edge',21);
electromagneticBC(pdem,'Voltage',0,'edge',27);
electromagneticBC(pdem,'Voltage',0,'edge',25);
electromagneticBC(pdem,'Voltage',0,'edge',31);
electromagneticBC(pdem,'Voltage',0,'edge',29);
electromagneticBC(pdem,'Voltage',0,'edge',35);
electromagneticBC(pdem,'Voltage',0,'edge',33);
electromagneticBC(pdem,'Voltage',0,'edge',39);
electromagneticBC(pdem,'Voltage',0,'edge',37);
electromagneticBC(pdem,'Voltage',0,'edge',43);
electromagneticBC(pdem,'Voltage',0,'edge',41);
% Negative up-strips
electromagneticBC(pdem,'Voltage',0,'edge',47);
electromagneticBC(pdem,'Voltage',0,'edge',45);
electromagneticBC(pdem,'Voltage',0,'edge',51);
electromagneticBC(pdem,'Voltage',0,'edge',49);
electromagneticBC(pdem,'Voltage',0,'edge',55);
electromagneticBC(pdem,'Voltage',0,'edge',53);
electromagneticBC(pdem,'Voltage',0,'edge',59);
electromagneticBC(pdem,'Voltage',0,'edge',57);
electromagneticBC(pdem,'Voltage',0,'edge',63);
electromagneticBC(pdem,'Voltage',0,'edge',61);
electromagneticBC(pdem,'Voltage',0,'edge',80);
electromagneticBC(pdem,'Voltage',0,'edge',65);
% Top all up-strips
electromagneticBC(pdem,'Voltage',0,'edge',14);
electromagneticBC(pdem,'Voltage',0,'edge',15);
electromagneticBC(pdem,'Voltage',0,'edge',16);
electromagneticBC(pdem,'Voltage',0,'edge',17);
electromagneticBC(pdem,'Voltage',0,'edge',18);
electromagneticBC(pdem,'Voltage',0,'edge',19);

electromagneticBC(pdem,'Voltage',BiasW,'edge',13);

electromagneticBC(pdem,'Voltage',0,'edge',12);
electromagneticBC(pdem,'Voltage',0,'edge',11);
electromagneticBC(pdem,'Voltage',0,'edge',10);
electromagneticBC(pdem,'Voltage',0,'edge',9);
electromagneticBC(pdem,'Voltage',0,'edge',8);
electromagneticBC(pdem,'Voltage',0,'edge',7);
% Bottom all up-strips
electromagneticBC(pdem,'Voltage',0,'edge',96);
electromagneticBC(pdem,'Voltage',0,'edge',98);
electromagneticBC(pdem,'Voltage',0,'edge',100);
electromagneticBC(pdem,'Voltage',0,'edge',102);
electromagneticBC(pdem,'Voltage',0,'edge',104);
electromagneticBC(pdem,'Voltage',0,'edge',106);

electromagneticBC(pdem,'Voltage',BiasW,'edge',94);

electromagneticBC(pdem,'Voltage',0,'edge',92);
electromagneticBC(pdem,'Voltage',0,'edge',90);
electromagneticBC(pdem,'Voltage',0,'edge',88);
electromagneticBC(pdem,'Voltage',0,'edge',86);
electromagneticBC(pdem,'Voltage',0,'edge',84);
electromagneticBC(pdem,'Voltage',0,'edge',82);
% Top edge
electromagneticBC(pdem,'Voltage',0,'edge',1);
% Right edge air
electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',137);
% Right edge sensor
electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',136);
% Right edge air
electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',135);
% Bottom edge
electromagneticBC(pdem,'Voltage',0,'edge',2);
% Left edge air
electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',140);
% Left edge sensor
electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',139);
% Left edge air
electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',138);
% Central down-strip
electromagneticBC(pdem,'Voltage',BiasV,'edge',3);
electromagneticBC(pdem,'Voltage',BiasV,'edge',5);
% Positive down-strips
electromagneticBC(pdem,'Voltage',BiasV,'edge',22);
electromagneticBC(pdem,'Voltage',BiasV,'edge',20);
electromagneticBC(pdem,'Voltage',BiasV,'edge',26);
electromagneticBC(pdem,'Voltage',BiasV,'edge',24);
electromagneticBC(pdem,'Voltage',BiasV,'edge',30);
electromagneticBC(pdem,'Voltage',BiasV,'edge',28);
electromagneticBC(pdem,'Voltage',BiasV,'edge',34);
electromagneticBC(pdem,'Voltage',BiasV,'edge',32);
electromagneticBC(pdem,'Voltage',BiasV,'edge',38);
electromagneticBC(pdem,'Voltage',BiasV,'edge',36);
electromagneticBC(pdem,'Voltage',BiasV,'edge',42);
electromagneticBC(pdem,'Voltage',BiasV,'edge',40);
% Negative down-strips
electromagneticBC(pdem,'Voltage',BiasV,'edge',46);
electromagneticBC(pdem,'Voltage',BiasV,'edge',44);
electromagneticBC(pdem,'Voltage',BiasV,'edge',50);
electromagneticBC(pdem,'Voltage',BiasV,'edge',48);
electromagneticBC(pdem,'Voltage',BiasV,'edge',54);
electromagneticBC(pdem,'Voltage',BiasV,'edge',52);
electromagneticBC(pdem,'Voltage',BiasV,'edge',58);
electromagneticBC(pdem,'Voltage',BiasV,'edge',56);
electromagneticBC(pdem,'Voltage',BiasV,'edge',62);
electromagneticBC(pdem,'Voltage',BiasV,'edge',60);
electromagneticBC(pdem,'Voltage',BiasV,'edge',79);
electromagneticBC(pdem,'Voltage',BiasV,'edge',64);
% Top all down-strips
electromagneticBC(pdem,'Voltage',BiasV,'edge',123);
electromagneticBC(pdem,'Voltage',BiasV,'edge',125);
electromagneticBC(pdem,'Voltage',BiasV,'edge',127);
electromagneticBC(pdem,'Voltage',BiasV,'edge',129);
electromagneticBC(pdem,'Voltage',BiasV,'edge',131);
electromagneticBC(pdem,'Voltage',BiasV,'edge',133);

electromagneticBC(pdem,'Voltage',BiasV,'edge',121);

electromagneticBC(pdem,'Voltage',BiasV,'edge',119);
electromagneticBC(pdem,'Voltage',BiasV,'edge',117);
electromagneticBC(pdem,'Voltage',BiasV,'edge',115);
electromagneticBC(pdem,'Voltage',BiasV,'edge',113);
electromagneticBC(pdem,'Voltage',BiasV,'edge',111);
electromagneticBC(pdem,'Voltage',BiasV,'edge',109);
% Bottom all down-strips
electromagneticBC(pdem,'Voltage',BiasV,'edge',73);
electromagneticBC(pdem,'Voltage',BiasV,'edge',74);
electromagneticBC(pdem,'Voltage',BiasV,'edge',75);
electromagneticBC(pdem,'Voltage',BiasV,'edge',76);
electromagneticBC(pdem,'Voltage',BiasV,'edge',77);
electromagneticBC(pdem,'Voltage',BiasV,'edge',78);

electromagneticBC(pdem,'Voltage',BiasV,'edge',72);

electromagneticBC(pdem,'Voltage',BiasV,'edge',71);
electromagneticBC(pdem,'Voltage',BiasV,'edge',70);
electromagneticBC(pdem,'Voltage',BiasV,'edge',69);
electromagneticBC(pdem,'Voltage',BiasV,'edge',68);
electromagneticBC(pdem,'Voltage',BiasV,'edge',67);
electromagneticBC(pdem,'Voltage',BiasV,'edge',66);


%%%%%%%%%%%%%%%%%
% Generate mesh %
%%%%%%%%%%%%%%%%%
msh = generateMesh(pdem,'Hmax',MeshMax,'GeometricOrder','quadratic');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%
pdem.VacuumPermittivity = eps0;
electromagneticSource(pdem,'face',1,'ChargeDensity',0);               % Air
electromagneticProperties(pdem,'RelativePermittivity',1,'face',1);    % Air
electromagneticSource(pdem,'face',3,'ChargeDensity',rho);             % Sensor
electromagneticProperties(pdem,'RelativePermittivity',epsR,'face',3); % Sensor
electromagneticSource(pdem,'face',2,'ChargeDensity',0);               % Air
electromagneticProperties(pdem,'RelativePermittivity',1,'face',2);    % Air
Potential = solve(pdem);


fprintf('CPU time --> %.2f [min]\n\n',(cputime-TStart)/60);
end
