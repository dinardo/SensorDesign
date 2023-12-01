%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation to compute the potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk   = Bulk thickness [um]
% PitchX = Pitch along X [um]
% BiasB  = Sensor backplane voltage [V] [0 Weighting; -V All]
% BiasW  = Sensor central strip voltage [V] [1 Weighting; 0 All]
% epsR   = Relative permittivity
% rho    = Charge density in the bulk [(Coulomb/um^3) / eps0 [F/um]]

function [pdem,Potential,DecomposedGeom,BulkStart,BulkStop] = StripPlanarMulti_SolvePoisson2D(...
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
NStrips      = 13;        % Total number of strips 
NStripLayers = 4;         % Number of strip layers


%%%%%%%%%%%%%%%%%%%%
% Create PDE model %
%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m solving Poisson equation in 2D to calculate the potential @@@\n');
pdem = createpde(1);


%%%%%%%%%%%%%%%%%%%%%%
% Create 2D geometry %
%%%%%%%%%%%%%%%%%%%%%%
gd = zeros(10,NStripLayers*NStrips+2);
sf = 'RWV - (';
ns = zeros(3,NStripLayers*NStrips+2);

for l = 1:NStripLayers
    % Central strips
    indx = 1+(l-1)*NStrips;
    gd(:,indx) = [ 3 4 -MetalWidth/2 MetalWidth/2 MetalWidth/2 -MetalWidth/2 ...
            BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers BulkStart+l*Bulk/NStripLayers ]';
    sf = strcat(sf,sprintf('R%02d+',indx));
    ns(:,indx) = sprintf('R%02d',indx);
    % Positive strips
    for s = 1:(NStrips-1)/2
        indx = s+1+(l-1)*NStrips;
        gd(:,indx) = [ 3 4 s*PitchX-MetalWidth/2 s*PitchX+MetalWidth/2 s*PitchX+MetalWidth/2 s*PitchX-MetalWidth/2 ...
            BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers BulkStart+l*Bulk/NStripLayers ]';
        sf = strcat(sf,sprintf('R%02d+',indx));
        ns(:,indx) = sprintf('R%02d',indx);
    end
    % Negative strips
    for s = 1:(NStrips-1)/2
        indx = s+1+(NStrips-1)/2+(l-1)*NStrips;
        gd(:,s+1+(NStrips-1)/2+(l-1)*NStrips) = [ 3 4 -(s*PitchX-MetalWidth/2) -(s*PitchX+MetalWidth/2) -(s*PitchX+MetalWidth/2) -(s*PitchX-MetalWidth/2) ...
            BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers+MetalThick BulkStart+l*Bulk/NStripLayers BulkStart+l*Bulk/NStripLayers ]';
        sf = strcat(sf,sprintf('R%02d+',indx));
        ns(:,indx) = sprintf('R%02d',indx);
    end
end
% Sensor volume RSV
gd(:,NStrips*NStripLayers+1) = [ 3 4 -((NStrips-1)/2*PitchX+PitchX/2) ((NStrips-1)/2*PitchX+PitchX/2) ((NStrips-1)/2*PitchX+PitchX/2) -((NStrips-1)/2*PitchX+PitchX/2) ...
    BulkStop BulkStop BulkStart BulkStart ]';
% Whole volume RWV
gd(:,NStrips*NStripLayers+2) = [ 3 4 -((NStrips-1)/2*PitchX+PitchX/2) ((NStrips-1)/2*PitchX+PitchX/2) ((NStrips-1)/2*PitchX+PitchX/2) -((NStrips-1)/2*PitchX+PitchX/2) ...
    Bulk*VolumeHeight Bulk*VolumeHeight 0 0 ]';

sf = strcat(sf(1:end-1),') + RSV');
ns(:,NStrips*NStripLayers+1) = 'RSV';
ns(:,NStrips*NStripLayers+2) = 'RWV';
ns = char(ns);

DecomposedGeom = decsg(gd,sf,ns);
geometryFromEdges(pdem,DecomposedGeom);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply boundary conditions (only on conductors) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All strips
applyBoundaryCondition(pdem,'neumann','edge',setdiff(1:pdem.Geometry.NumEdges,...
    linspace(28,28+NStrips*2,NStrips+1)),'q',0,'g',0);
% Top strips
applyBoundaryCondition(pdem,'dirichlet','edge',15:15+NStrips-1,'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',linspace(29,29+(NStrips-1)*2,NStrips),'h',1,'r',0);
applyBoundaryCondition(pdem,'dirichlet','edge',[213,195,165,148,144,140,136,132,...
    128,124,120,116,64,60,72,68,80,76,88,84,98,92,104,100,112,108],'h',1,'r',0);
% Central strip
applyBoundaryCondition(pdem,'dirichlet','edge',21,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'dirichlet','edge',41,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'dirichlet','edge',60,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'dirichlet','edge',64,'h',1,'r',BiasW);
% Top edge
applyBoundaryCondition(pdem,'dirichlet','edge',209,'h',1,'r',0);
% Right edge sensor
applyBoundaryCondition(pdem,'neumann','edge',55,'q',0,'g',0);
% Right edge air
applyBoundaryCondition(pdem,'neumann','edge',56,'q',0,'g',0);
% Left edge sensor
applyBoundaryCondition(pdem,'neumann','edge',227,'q',0,'g',0);
% Left edge air
applyBoundaryCondition(pdem,'neumann','edge',228,'q',0,'g',0);
% Bottom edge
applyBoundaryCondition(pdem,'dirichlet','edge',1,'h',1,'r',BiasB);


%%%%%%%%%%%%%%%%%
% Generate mesh %
%%%%%%%%%%%%%%%%%
msh = generateMesh(pdem,'Hmax',MeshMax,'GeometricOrder','quadratic');


% @TMP@
myXlim = [-PitchX * NStrips/2,+PitchX * NStrips/2];
myYlim = [0,Bulk * VolumeHeight];
figure(1);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%
specifyCoefficients(pdem,'m',0,'d',0,'c',1,   'a',0,'f',0,  'face',2); % Air
specifyCoefficients(pdem,'m',0,'d',0,'c',epsR,'a',0,'f',rho,'face',1); % Sensor
Potential = solvepde(pdem);


fprintf('CPU time --> %.2f [min]\n\n',(cputime-TStart)/60);
end
