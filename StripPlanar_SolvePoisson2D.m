%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation to compute the potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk   = Bulk thickness [um]
% PitchX = Pitch along X [um]
% BiasV  = Sensor backplane voltage [V] == 0 ? compute weighting field
% epsR   = Relative permittivity
% rho    = Charge density in the bulk [(Coulomb/um^3)]

function [pdem,Potential,DecomposedGeom,BulkStart,BulkStop,VolumeHeight] =...
    StripPlanar_SolvePoisson2D(Bulk,PitchX,BiasV,epsR,rho)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
NStripLayers = 1;         % Number of strip layers
eps0         = 8.85e-18;  % Vacuum permittivity [F/um]
VolumeHeight = 3;         % Volume height [units of bulk thickness]
MeshMax      = 5;         % Maximum mesh edge length [um]
MetalThick   = 5;         % Metalization thickness [um]
MetalWidth   = PitchX-20; % Metalization width [um]
BulkStart    = 0;         % Bulk start coordinate [um]
BulkStop     = Bulk+(NStripLayers-1)*MetalThick; % Bulk stop coordinate [um]
NStrips      = 13;        % Total number of strips 
BiasW        = 0;         % Bias to compute weighting potential
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
gd         = zeros(10,NStripLayers*NStrips+2);
bulkStrips = ' (RSV - (';
sf         = 'RWV - (';
ns         = zeros(3,NStripLayers*NStrips+2);

for l = 1:NStripLayers
    % Central strips
    indx = 1+(l-1)*NStrips;
    gd(:,indx) = [ 3 4 -MetalWidth/2 MetalWidth/2 MetalWidth/2 -MetalWidth/2 ...
            BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers BulkStart+l*Bulk/NStripLayers ]';
    sf = strcat(sf,sprintf('R%02d+',indx));
    ns(:,indx) = sprintf('R%02d',indx);
    bulkStrips = strcat(bulkStrips,sprintf('R%02d+',indx));
    % Positive strips
    for s = 1:(NStrips-1)/2
        indx = s+1+(l-1)*NStrips;
        gd(:,indx) = [ 3 4 s*PitchX-MetalWidth/2 s*PitchX+MetalWidth/2 s*PitchX+MetalWidth/2 s*PitchX-MetalWidth/2 ...
            BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers BulkStart+l*Bulk/NStripLayers ]';
        sf = strcat(sf,sprintf('R%02d+',indx));
        ns(:,indx) = sprintf('R%02d',indx);
        bulkStrips = strcat(bulkStrips,sprintf('R%02d+',indx));
    end
    % Negative strips
    for s = 1:(NStrips-1)/2
        indx = s+1+(NStrips-1)/2+(l-1)*NStrips;
        gd(:,s+1+(NStrips-1)/2+(l-1)*NStrips) = [ 3 4 -(s*PitchX-MetalWidth/2) -(s*PitchX+MetalWidth/2) -(s*PitchX+MetalWidth/2) -(s*PitchX-MetalWidth/2) ...
            BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers BulkStart+l*Bulk/NStripLayers ]';
        sf = strcat(sf,sprintf('R%02d+',indx));
        ns(:,indx) = sprintf('R%02d',indx);
        bulkStrips = strcat(bulkStrips,sprintf('R%02d+',indx));
    end
end
% Sensor volume RSV
gd(:,NStrips*NStripLayers+1) = [ 3 4 -((NStrips-1)/2*PitchX+PitchX/2) ((NStrips-1)/2*PitchX+PitchX/2) ((NStrips-1)/2*PitchX+PitchX/2) -((NStrips-1)/2*PitchX+PitchX/2) ...
    BulkStop BulkStop BulkStart BulkStart ]';
% Whole volume RWV
gd(:,NStrips*NStripLayers+2) = [ 3 4 -((NStrips-1)/2*PitchX+PitchX/2) ((NStrips-1)/2*PitchX+PitchX/2) ((NStrips-1)/2*PitchX+PitchX/2) -((NStrips-1)/2*PitchX+PitchX/2) ...
    Bulk*VolumeHeight Bulk*VolumeHeight BulkStart BulkStart ]';

bulkStrips = strcat(bulkStrips(1:end-1),'))');
sf = strcat(sf(1:end-1),strcat(') +',bulkStrips));
ns(:,NStrips*NStripLayers+1) = 'RSV';
ns(:,NStrips*NStripLayers+2) = 'RWV';
ns = char(ns);

DecomposedGeom = decsg(gd,sf,ns);
geometryFromEdges(pdem,DecomposedGeom);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply boundary conditions (only on conductors) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if NStripLayers == 1
    % Central strip
    electromagneticBC(pdem,'Voltage',BiasW,'edge',1:2);
    % Positive strips
    electromagneticBC(pdem,'Voltage',0,'edge',3:14);
    % Negative strips
    electromagneticBC(pdem,'Voltage',0,'edge',15:26);
    % Top all strips
    electromagneticBC(pdem,'Voltage',0,'edge',28:33);
    electromagneticBC(pdem,'Voltage',BiasW,'edge',34);
    electromagneticBC(pdem,'Voltage',0,'edge',35:40);
    % Bottom all strips
    electromagneticBC(pdem,'Voltage',0,'edge',linspace(42,52,6));
    electromagneticBC(pdem,'Voltage',BiasW,'edge',54);
    electromagneticBC(pdem,'Voltage',0,'edge',linspace(56,66,6));
    % Top edge
    electromagneticBC(pdem,'Voltage',0,'edge',27);
    % Right edge sensor
    electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',68);
    % Right edge air
    electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',69);
    % Bottom edge
    electromagneticBC(pdem,'Voltage',BiasV,'edge',70);
    % Left edge sensor
    electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',71);
    % Left edge air
    electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',72);
else
    % All internal strips
    top    = [149 166 214];
    bottom = [179 196 2];
    for i = 0:NStripLayers-2
        electromagneticBC(pdem,'Voltage',BiasV/NStripLayers*(NStripLayers-i-1),'edge',top(i+1):top(i+1)+NStrips-1);
        electromagneticBC(pdem,'Voltage',BiasV/NStripLayers*(NStripLayers-i-1),'edge',bottom(i+1):bottom(i+1)+NStrips-1);
        electromagneticBC(pdem,'Voltage',BiasV/NStripLayers*(NStripLayers-i-1),'edge',linspace(61+i,109+i,(NStrips-1)/2+1));
        electromagneticBC(pdem,'Voltage',BiasV/NStripLayers*(NStripLayers-i-1),'edge',linspace(57+i,105+i,(NStrips-1)/2+1));
        electromagneticBC(pdem,'Voltage',BiasV/NStripLayers*(NStripLayers-i-1),'edge',linspace(113+i,145+i,(NStrips-1)/2-1));
        electromagneticBC(pdem,'Voltage',BiasV/NStripLayers*(NStripLayers-i-1),'edge',linspace(117+i,141+i,(NStrips-1)/2-2));
        electromagneticBC(pdem,'Voltage',BiasV/NStripLayers*(NStripLayers-i-1),'edge',192+i);
        electromagneticBC(pdem,'Voltage',BiasV/NStripLayers*(NStripLayers-i-1),'edge',210+i);
        electromagneticBC(pdem,'Voltage',BiasV/NStripLayers*(NStripLayers-i-1),'edge',162+i);
    end
    % Top strips
    electromagneticBC(pdem,'Voltage',0,'edge',15:15+NStrips-1);
    electromagneticBC(pdem,'Voltage',0,'edge',linspace(29,29+(NStrips-1)*2,NStrips));
    electromagneticBC(pdem,'Voltage',0,'edge',[213,195,165,148,144,140,136,132,...
        128,124,120,116,72,68,80,76,88,84,96,92,104,100,112,108]);
    % Central top strip
    electromagneticBC(pdem,'Voltage',BiasW,'edge',21);
    electromagneticBC(pdem,'Voltage',BiasW,'edge',41);
    electromagneticBC(pdem,'Voltage',BiasW,'edge',60);
    electromagneticBC(pdem,'Voltage',BiasW,'edge',64);
    % Top edge
    electromagneticBC(pdem,'Voltage',0,'edge',209);
    % Right edge sensor
    electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',55);
    % Right edge air
    electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',56);
    % Left edge sensor
    electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',227);
    % Left edge air
    electromagneticBC(pdem,'SurfaceCurrentDensity',0,'edge',228);
    % Bottom edge
    electromagneticBC(pdem,'Voltage',BiasV,'edge',1);
end


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
electromagneticSource(pdem,'face',2,'ChargeDensity',rho);             % Sensor
electromagneticProperties(pdem,'RelativePermittivity',epsR,'face',2); % Sensor
Potential = solve(pdem);


fprintf('CPU time --> %.2f [min]\n\n',(cputime-TStart)/60);
end
