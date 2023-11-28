%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation to compute the potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk    = Bulk thickness [um]
% PitchX  = Pitch along X [um]
% PitchY  = Pitch along Y [um]
% BiasB   = Sensor backplane voltage [V] [0 Weighting; -V All]
% BiasW   = Sensor central pixel voltage [V] [1 Weighting; 0 All]
% epsR    = Relative permittivity
% rho     = Charge density in the bulk [(Coulomb/um^3) / eps0 [F/um]]
% XQ      = Coordinate for potential query along z [um]
% YQ      = Coordinate for potential query along z [um]
% ItFigIn = Figure iterator input

function [Potential, Sq, zq, ItFigOut] = PixelPlanar_SolvePoisson3D(...
    Bulk,PitchX,PitchY,BiasB,BiasW,epsR,rho,XQ,YQ,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
VolumeHeight = 2;         % Volume height [units of bulk thickness]
MetalThick   = 5;         % Metalization thickness [um]
MetalWidthX  = PitchX-20; % Metalization width along X [um]
MetalWidthY  = PitchY-20; % Metalization width along Y [um]
ReSampleFine = 1;         % Used in order to make nice plots [um]
StepMeshHol  = 1;         % Step to build mesh hollow volume [um]
StepMeshVol  = 4;         % Step to build mesh whole volume [um]
StepSlices   = Bulk/10;   % Step to build slices along z [um]


%%%%%%%%%%%%%%%%%%%%
% Create PDE model %
%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m solving Poisson equation in 3D to calculate the potential @@@\n');
pdem = createpde(1);


%%%%%%%%%%%%%%%%%%%%%%
% Create 3D geometry %
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% Row 0 %
%%%%%%%%%
[x00,y00,z00] = meshgrid(-MetalWidthX/2:StepMeshHol:MetalWidthX/2,...
    -MetalWidthY/2:StepMeshHol:MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

x00 = x00(:);
y00 = y00(:);
z00 = z00(:);

frontier = (x00 == -MetalWidthX/2 | x00 == MetalWidthX/2 | ...
            y00 == -MetalWidthY/2 | y00 == MetalWidthY/2 | ...
            z00 == Bulk           | z00 == Bulk+MetalThick);

x00(frontier) = [];
y00(frontier) = [];
z00(frontier) = [];

[xp10,yp10,zp10] = meshgrid(1*PitchX-MetalWidthX/2:StepMeshHol:1*PitchX+MetalWidthX/2,...
    -MetalWidthY/2:StepMeshHol:MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xp10 = xp10(:);
yp10 = yp10(:);
zp10 = zp10(:);

frontier = (xp10 == 1*PitchX-MetalWidthX/2 | xp10 == 1*PitchX+MetalWidthX/2 | ...
            yp10 == -MetalWidthY/2         | yp10 == MetalWidthY/2          | ...
            zp10 == Bulk                   | zp10 == Bulk+MetalThick);

xp10(frontier) = [];
yp10(frontier) = [];
zp10(frontier) = [];

[xp20,yp20,zp20] = meshgrid(2*PitchX-MetalWidthX/2:StepMeshHol:2*PitchX+MetalWidthX/2,...
    -MetalWidthY/2:StepMeshHol:MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xp20 = xp20(:);
yp20 = yp20(:);
zp20 = zp20(:);

frontier = (xp20 == 2*PitchX-MetalWidthX/2 | xp20 == 2*PitchX+MetalWidthX/2 | ...
            yp20 == -MetalWidthY/2         | yp20 == MetalWidthY/2          | ...
            zp20 == Bulk                   | zp20 == Bulk+MetalThick);

xp20(frontier) = [];
yp20(frontier) = [];
zp20(frontier) = [];

[xm10,ym10,zm10] = meshgrid(-1*PitchX-MetalWidthX/2:StepMeshHol:-1*PitchX+MetalWidthX/2,...
    -MetalWidthY/2:StepMeshHol:MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xm10 = xm10(:);
ym10 = ym10(:);
zm10 = zm10(:);

frontier = (xm10 == -1*PitchX-MetalWidthX/2 | xm10 == -1*PitchX+MetalWidthX/2 | ...
            ym10 == -MetalWidthY/2          | ym10 == MetalWidthY/2           | ...
            zm10 == Bulk                    | zm10 == Bulk+MetalThick);

xm10(frontier) = [];
ym10(frontier) = [];
zm10(frontier) = [];

[xm20,ym20,zm20] = meshgrid(-2*PitchX-MetalWidthX/2:StepMeshHol:-2*PitchX+MetalWidthX/2,...
    -MetalWidthY/2:StepMeshHol:MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xm20 = xm20(:);
ym20 = ym20(:);
zm20 = zm20(:);

frontier = (xm20 == -2*PitchX-MetalWidthX/2 | xm20 == -2*PitchX+MetalWidthX/2 | ...
            ym20 == -MetalWidthY/2          | ym20 == MetalWidthY/2           | ...
            zm20 == Bulk                    | zm20 == Bulk+MetalThick);

xm20(frontier) = [];
ym20(frontier) = [];
zm20(frontier) = [];

%%%%%%%%%%
% Row +1 %
%%%%%%%%%%
[x0p1,y0p1,z0p1] = meshgrid(-MetalWidthX/2:StepMeshHol:MetalWidthX/2,...
    1*PitchY-MetalWidthY/2:StepMeshHol:1*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

x0p1 = x0p1(:);
y0p1 = y0p1(:);
z0p1 = z0p1(:);

frontier = (x0p1 == -MetalWidthX/2         | x0p1 == MetalWidthX/2          | ...
            y0p1 == 1*PitchY-MetalWidthY/2 | y0p1 == 1*PitchY+MetalWidthY/2 | ...
            z0p1 == Bulk                   | z0p1 == Bulk+MetalThick);

x0p1(frontier) = [];
y0p1(frontier) = [];
z0p1(frontier) = [];

[xp1p1,yp1p1,zp1p1] = meshgrid(1*PitchX-MetalWidthX/2:StepMeshHol:1*PitchX+MetalWidthX/2,...
    1*PitchY-MetalWidthY/2:StepMeshHol:1*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xp1p1 = xp1p1(:);
yp1p1 = yp1p1(:);
zp1p1 = zp1p1(:);

frontier = (xp1p1 == 1*PitchX-MetalWidthX/2 | xp1p1 == 1*PitchX+MetalWidthX/2 | ...
            yp1p1 == 1*PitchY-MetalWidthY/2 | yp1p1 == 1*PitchY+MetalWidthY/2 | ...
            zp1p1 == Bulk                   | zp1p1 == Bulk+MetalThick);

xp1p1(frontier) = [];
yp1p1(frontier) = [];
zp1p1(frontier) = [];

[xp2p1,yp2p1,zp2p1] = meshgrid(2*PitchX-MetalWidthX/2:StepMeshHol:2*PitchX+MetalWidthX/2,...
    1*PitchY-MetalWidthY/2:StepMeshHol:1*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xp2p1 = xp2p1(:);
yp2p1 = yp2p1(:);
zp2p1 = zp2p1(:);

frontier = (xp2p1 == 2*PitchX-MetalWidthX/2 | xp2p1 == 2*PitchX+MetalWidthX/2 | ...
            yp2p1 == 1*PitchY-MetalWidthY/2 | yp2p1 == 1*PitchY+MetalWidthY/2 | ...
            zp2p1 == Bulk                   | zp2p1 == Bulk+MetalThick);

xp2p1(frontier) = [];
yp2p1(frontier) = [];
zp2p1(frontier) = [];

[xm1p1,ym1p1,zm1p1] = meshgrid(-1*PitchX-MetalWidthX/2:StepMeshHol:-1*PitchX+MetalWidthX/2,...
    1*PitchY-MetalWidthY/2:StepMeshHol:1*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xm1p1 = xm1p1(:);
ym1p1 = ym1p1(:);
zm1p1 = zm1p1(:);

frontier = (xm1p1 == -1*PitchX-MetalWidthX/2 | xm1p1 == -1*PitchX+MetalWidthX/2 | ...
            ym1p1 == 1*PitchY-MetalWidthY/2  | ym1p1 == 1*PitchY+MetalWidthY/2  | ...
            zm1p1 == Bulk                    | zm1p1 == Bulk+MetalThick);

xm1p1(frontier) = [];
ym1p1(frontier) = [];
zm1p1(frontier) = [];

[xm2p1,ym2p1,zm2p1] = meshgrid(-2*PitchX-MetalWidthX/2:StepMeshHol:-2*PitchX+MetalWidthX/2,...
    1*PitchY-MetalWidthY/2:StepMeshHol:1*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xm2p1 = xm2p1(:);
ym2p1 = ym2p1(:);
zm2p1 = zm2p1(:);

frontier = (xm2p1 == -2*PitchX-MetalWidthX/2 | xm2p1 == -2*PitchX+MetalWidthX/2 | ...
            ym2p1 == 1*PitchY-MetalWidthY/2  | ym2p1 == 1*PitchY+MetalWidthY/2  | ...
            zm2p1 == Bulk                    | zm2p1 == Bulk+MetalThick);

xm2p1(frontier) = [];
ym2p1(frontier) = [];
zm2p1(frontier) = [];

%%%%%%%%%%
% Row +2 %
%%%%%%%%%%
[x0p2,y0p2,z0p2] = meshgrid(-MetalWidthX/2:StepMeshHol:MetalWidthX/2,...
    2*PitchY-MetalWidthY/2:StepMeshHol:2*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

x0p2 = x0p2(:);
y0p2 = y0p2(:);
z0p2 = z0p2(:);

frontier = (x0p2 == -MetalWidthX/2         | x0p2 == MetalWidthX/2          | ...
            y0p2 == 2*PitchY-MetalWidthY/2 | y0p2 == 2*PitchY+MetalWidthY/2 | ...
            z0p2 == Bulk                   | z0p2 == Bulk+MetalThick);

x0p2(frontier) = [];
y0p2(frontier) = [];
z0p2(frontier) = [];

[xp1p2,yp1p2,zp1p2] = meshgrid(1*PitchX-MetalWidthX/2:StepMeshHol:1*PitchX+MetalWidthX/2,...
    2*PitchY-MetalWidthY/2:StepMeshHol:2*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xp1p2 = xp1p2(:);
yp1p2 = yp1p2(:);
zp1p2 = zp1p2(:);

frontier = (xp1p2 == 1*PitchX-MetalWidthX/2 | xp1p2 == 1*PitchX+MetalWidthX/2 | ...
            yp1p2 == 2*PitchY-MetalWidthY/2 | yp1p2 == 2*PitchY+MetalWidthY/2 | ...
            zp1p2 == Bulk                   | zp1p2 == Bulk+MetalThick);

xp1p2(frontier) = [];
yp1p2(frontier) = [];
zp1p2(frontier) = [];

[xp2p2,yp2p2,zp2p2] = meshgrid(2*PitchX-MetalWidthX/2:StepMeshHol:2*PitchX+MetalWidthX/2,...
    2*PitchY-MetalWidthY/2:StepMeshHol:2*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xp2p2 = xp2p2(:);
yp2p2 = yp2p2(:);
zp2p2 = zp2p2(:);

frontier = (xp2p2 == 2*PitchX-MetalWidthX/2 | xp2p2 == 2*PitchX+MetalWidthX/2 | ...
            yp2p2 == 2*PitchY-MetalWidthY/2 | yp2p2 == 2*PitchY+MetalWidthY/2 | ...
            zp2p2 == Bulk                   | zp2p2 == Bulk+MetalThick);

xp2p2(frontier) = [];
yp2p2(frontier) = [];
zp2p2(frontier) = [];

[xm1p2,ym1p2,zm1p2] = meshgrid(-1*PitchX-MetalWidthX/2:StepMeshHol:-1*PitchX+MetalWidthX/2,...
    2*PitchY-MetalWidthY/2:StepMeshHol:2*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xm1p2 = xm1p2(:);
ym1p2 = ym1p2(:);
zm1p2 = zm1p2(:);

frontier = (xm1p2 == -1*PitchX-MetalWidthX/2 | xm1p2 == -1*PitchX+MetalWidthX/2 | ...
            ym1p2 == 2*PitchY-MetalWidthY/2  | ym1p2 == 2*PitchY+MetalWidthY/2  | ...
            zm1p2 == Bulk                    | zm1p2 == Bulk+MetalThick);

xm1p2(frontier) = [];
ym1p2(frontier) = [];
zm1p2(frontier) = [];

[xm2p2,ym2p2,zm2p2] = meshgrid(-2*PitchX-MetalWidthX/2:StepMeshHol:-2*PitchX+MetalWidthX/2,...
    2*PitchY-MetalWidthY/2:StepMeshHol:2*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xm2p2 = xm2p2(:);
ym2p2 = ym2p2(:);
zm2p2 = zm2p2(:);

frontier = (xm2p2 == -2*PitchX-MetalWidthX/2 | xm2p2 == -2*PitchX+MetalWidthX/2 | ...
            ym2p2 == 2*PitchY-MetalWidthY/2  | ym2p2 == 2*PitchY+MetalWidthY/2  | ...
            zm2p2 == Bulk                    | zm2p2 == Bulk+MetalThick);

xm2p2(frontier) = [];
ym2p2(frontier) = [];
zm2p2(frontier) = [];

%%%%%%%%%%
% Row -1 %
%%%%%%%%%%
[x0m1,y0m1,z0m1] = meshgrid(-MetalWidthX/2:StepMeshHol:MetalWidthX/2,...
    -1*PitchY-MetalWidthY/2:StepMeshHol:-1*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

x0m1 = x0m1(:);
y0m1 = y0m1(:);
z0m1 = z0m1(:);

frontier = (x0m1 == -MetalWidthX/2          | x0m1 == MetalWidthX/2           | ...
            y0m1 == -1*PitchY-MetalWidthY/2 | y0m1 == -1*PitchY+MetalWidthY/2 | ...
            z0m1 == Bulk                    | z0m1 == Bulk+MetalThick);

x0m1(frontier) = [];
y0m1(frontier) = [];
z0m1(frontier) = [];

[xp1m1,yp1m1,zp1m1] = meshgrid(1*PitchX-MetalWidthX/2:StepMeshHol:1*PitchX+MetalWidthX/2,...
    -1*PitchY-MetalWidthY/2:StepMeshHol:-1*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xp1m1 = xp1m1(:);
yp1m1 = yp1m1(:);
zp1m1 = zp1m1(:);

frontier = (xp1m1 == 1*PitchX-MetalWidthX/2  | xp1m1 == 1*PitchX+MetalWidthX/2  | ...
            yp1m1 == -1*PitchY-MetalWidthY/2 | yp1m1 == -1*PitchY+MetalWidthY/2 | ...
            zp1m1 == Bulk                    | zp1m1 == Bulk+MetalThick);

xp1m1(frontier) = [];
yp1m1(frontier) = [];
zp1m1(frontier) = [];

[xp2m1,yp2m1,zp2m1] = meshgrid(2*PitchX-MetalWidthX/2:StepMeshHol:2*PitchX+MetalWidthX/2,...
    -1*PitchY-MetalWidthY/2:StepMeshHol:-1*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xp2m1 = xp2m1(:);
yp2m1 = yp2m1(:);
zp2m1 = zp2m1(:);

frontier = (xp2m1 == 2*PitchX-MetalWidthX/2  | xp2m1 == 2*PitchX+MetalWidthX/2  | ...
            yp2m1 == -1*PitchY-MetalWidthY/2 | yp2m1 == -1*PitchY+MetalWidthY/2 | ...
            zp2m1 == Bulk                    | zp2m1 == Bulk+MetalThick);

xp2m1(frontier) = [];
yp2m1(frontier) = [];
zp2m1(frontier) = [];

[xm1m1,ym1m1,zm1m1] = meshgrid(-1*PitchX-MetalWidthX/2:StepMeshHol:-1*PitchX+MetalWidthX/2,...
    -1*PitchY-MetalWidthY/2:StepMeshHol:-1*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xm1m1 = xm1m1(:);
ym1m1 = ym1m1(:);
zm1m1 = zm1m1(:);

frontier = (xm1m1 == -1*PitchX-MetalWidthX/2 | xm1m1 == -1*PitchX+MetalWidthX/2 | ...
            ym1m1 == -1*PitchY-MetalWidthY/2 | ym1m1 == -1*PitchY+MetalWidthY/2 | ...
            zm1m1 == Bulk                    | zm1m1 == Bulk+MetalThick);

xm1m1(frontier) = [];
ym1m1(frontier) = [];
zm1m1(frontier) = [];

[xm2m1,ym2m1,zm2m1] = meshgrid(-2*PitchX-MetalWidthX/2:StepMeshHol:-2*PitchX+MetalWidthX/2,...
    -1*PitchY-MetalWidthY/2:StepMeshHol:-1*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xm2m1 = xm2m1(:);
ym2m1 = ym2m1(:);
zm2m1 = zm2m1(:);

frontier = (xm2m1 == -2*PitchX-MetalWidthX/2 | xm2m1 == -2*PitchX+MetalWidthX/2 | ...
            ym2m1 == -1*PitchY-MetalWidthY/2 | ym2m1 == -1*PitchY+MetalWidthY/2 | ...
            zm2m1 == Bulk                    | zm2m1 == Bulk+MetalThick);

xm2m1(frontier) = [];
ym2m1(frontier) = [];
zm2m1(frontier) = [];

%%%%%%%%%%
% Row -2 %
%%%%%%%%%%
[x0m2,y0m2,z0m2] = meshgrid(-MetalWidthX/2:StepMeshHol:MetalWidthX/2,...
    -2*PitchY-MetalWidthY/2:StepMeshHol:-2*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

x0m2 = x0m2(:);
y0m2 = y0m2(:);
z0m2 = z0m2(:);

frontier = (x0m2 == -MetalWidthX/2          | x0m2 == MetalWidthX/2           | ...
            y0m2 == -2*PitchY-MetalWidthY/2 | y0m2 == -2*PitchY+MetalWidthY/2 | ...
            z0m2 == Bulk                    | z0m2 == Bulk+MetalThick);

x0m2(frontier) = [];
y0m2(frontier) = [];
z0m2(frontier) = [];

[xp1m2,yp1m2,zp1m2] = meshgrid(1*PitchX-MetalWidthX/2:StepMeshHol:1*PitchX+MetalWidthX/2,...
    -2*PitchY-MetalWidthY/2:StepMeshHol:-2*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xp1m2 = xp1m2(:);
yp1m2 = yp1m2(:);
zp1m2 = zp1m2(:);

frontier = (xp1m2 == 1*PitchX-MetalWidthX/2  | xp1m2 == 1*PitchX+MetalWidthX/2  | ...
            yp1m2 == -2*PitchY-MetalWidthY/2 | yp1m2 == -2*PitchY+MetalWidthY/2 | ...
            zp1m2 == Bulk                    | zp1m2 == Bulk+MetalThick);

xp1m2(frontier) = [];
yp1m2(frontier) = [];
zp1m2(frontier) = [];

[xp2m2,yp2m2,zp2m2] = meshgrid(2*PitchX-MetalWidthX/2:StepMeshHol:2*PitchX+MetalWidthX/2,...
    -2*PitchY-MetalWidthY/2:StepMeshHol:-2*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xp2m2 = xp2m2(:);
yp2m2 = yp2m2(:);
zp2m2 = zp2m2(:);

frontier = (xp2m2 == 2*PitchX-MetalWidthX/2  | xp2m2 == 2*PitchX+MetalWidthX/2  | ...
            yp2m2 == -2*PitchY-MetalWidthY/2 | yp2m2 == -2*PitchY+MetalWidthY/2 | ...
            zp2m2 == Bulk                    | zp2m2 == Bulk+MetalThick);

xp2m2(frontier) = [];
yp2m2(frontier) = [];
zp2m2(frontier) = [];

[xm1m2,ym1m2,zm1m2] = meshgrid(-1*PitchX-MetalWidthX/2:StepMeshHol:-1*PitchX+MetalWidthX/2,...
    -2*PitchY-MetalWidthY/2:StepMeshHol:-2*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xm1m2 = xm1m2(:);
ym1m2 = ym1m2(:);
zm1m2 = zm1m2(:);

frontier = (xm1m2 == -1*PitchX-MetalWidthX/2 | xm1m2 == -1*PitchX+MetalWidthX/2 | ...
            ym1m2 == -2*PitchY-MetalWidthY/2 | ym1m2 == -2*PitchY+MetalWidthY/2 | ...
            zm1m2 == Bulk                    | zm1m2 == Bulk+MetalThick);

xm1m2(frontier) = [];
ym1m2(frontier) = [];
zm1m2(frontier) = [];

[xm2m2,ym2m2,zm2m2] = meshgrid(-2*PitchX-MetalWidthX/2:StepMeshHol:-2*PitchX+MetalWidthX/2,...
    -2*PitchY-MetalWidthY/2:StepMeshHol:-2*PitchY+MetalWidthY/2,...
    Bulk:StepMeshHol:Bulk+MetalThick);

xm2m2 = xm2m2(:);
ym2m2 = ym2m2(:);
zm2m2 = zm2m2(:);

frontier = (xm2m2 == -2*PitchX-MetalWidthX/2 | xm2m2 == -2*PitchX+MetalWidthX/2 | ...
            ym2m2 == -2*PitchY-MetalWidthY/2 | ym2m2 == -2*PitchY+MetalWidthY/2 | ...
            zm2m2 == Bulk                    | zm2m2 == Bulk+MetalThick);

xm2m2(frontier) = [];
ym2m2(frontier) = [];
zm2m2(frontier) = [];

%%%%%%%%%%%%%%%%
% Whole volume %
%%%%%%%%%%%%%%%%
[x,y,z] = meshgrid(-2*PitchX-MetalWidthX:StepMeshVol:2*PitchX+MetalWidthX,...
    -2*PitchY-MetalWidthY:StepMeshVol:2*PitchY+MetalWidthY,...
    0:StepMeshVol:Bulk*VolumeHeight);

x = x(:);
y = y(:);
z = z(:);


%%%%%%%%%%%%%%%%%%%%%%
% Combine all pixels %
%%%%%%%%%%%%%%%%%%%%%%
xPx = [xm2m2' xm1m2' x0m2' xp1m2' xp2m2'...
       xm2m1' xm1m1' x0m1' xp1m1' xp2m1'...
       xm20'  xm10'  x00'  xp10'  xp20'...
       xm2p1' xm1p1' x0p1' xp1p1' xp2p1'...
       xm2p2' xm1p2' x0p2' xp1p2' xp2p2']';
yPx = [ym2m2' ym1m2' y0m2' yp1m2' yp2m2'...
       ym2m1' ym1m1' y0m1' yp1m1' yp2m1'...
       ym20'  ym10'  y00'  yp10'  yp20'...
       ym2p1' ym1p1' y0p1' yp1p1' yp2p1'...
       ym2p2' ym1p2' y0p2' yp1p2' yp2p2']';
zPx = [zm2m2' zm1m2' z0m2' zp1m2' zp2m2'...
       zm2m1' zm1m1' z0m1' zp1m1' zp2m1'...
       zm20'  zm10'  z00'  zp10'  zp20'...
       zm2p1' zm1p1' z0p1' zp1p1' zp2p1'...
       zm2p2' zm1p2' z0p2' zp1p2' zp2p2']';

shpPx = alphaShape(xPx,yPx,zPx);


%%%%%%%%%%%%%%%%%%%%%%%
% Create final volume %
%%%%%%%%%%%%%%%%%%%%%%%
inVol = inShape(shpPx,x,y,z);
x = x(~inVol);
y = y(~inVol);
z = z(~inVol);

shpAll = alphaShape(x,y,z);
[elements,nodes] = boundaryFacets(shpAll);
nodes = nodes';
elements = elements';


%%%%%%%%%%%%%%%%%%%%%%%%%
% Create final geometry %
%%%%%%%%%%%%%%%%%%%%%%%%%
geometryFromMesh(pdem,nodes,elements);


%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation for all faces
applyBoundaryCondition(pdem,'dirichlet','face',1:pdem.Geometry.NumFaces,'h',1,'r',0);
% Backplane
applyBoundaryCondition(pdem,'dirichlet','face',1,'h',1,'r',BiasB);
% Central pixel
applyBoundaryCondition(pdem,'dirichlet','face',22,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'dirichlet','face',35,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'dirichlet','face',44,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'dirichlet','face',80,'h',1,'r',BiasW);


%%%%%%%%%%%%%%%
% Create mesh %
%%%%%%%%%%%%%%%
msh = generateMesh(pdem,'GeometricOrder','quadratic');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%
specifyCoefficients(pdem,'m',0,'d',0,'c',epsR,'a',0,'f',rho);
Potential = solvepde(pdem);


%%%%%%%%%
% Plots %
%%%%%%%%%
figure(ItFigIn);
h = pdegplot(pdem,'FaceLabels','on');
h(1).FaceAlpha = 0.5;
title('Geometry');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Z [\mum]');

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
u = Potential.NodalSolution;
pdeplot3D(pdem,'colormapdata',u);
axis on;
title('Potential on conductors');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Z [\mum]');

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
colormap jet;
[x,y,z] = meshgrid(-2*PitchX-MetalWidthX/2:StepMeshVol:2*PitchX+MetalWidthX/2,...
    -2*PitchY-MetalWidthY/2:StepMeshVol:2*PitchY+MetalWidthY/2,...
    0:StepMeshVol:Bulk*3/2);
sl = interpolateSolution(Potential,x,y,z);
sl = reshape(sl,size(x));
contourslice(x,y,z,sl,[],[],0:StepSlices:Bulk);
title('Potential slices');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Z [\mum]');
grid on;
colorbar;

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
zq = 0:ReSampleFine:Bulk;
xq = XQ * ones(1,length(zq));
yq = YQ * ones(1,length(zq));
Sq = interpolateSolution(Potential,xq,yq,zq);
plot(zq,Sq);
title(sprintf('Potential along z at x = %.2f um y = %.2f um',XQ,YQ));
xlabel('Z [\mum]');
ylabel('Potential');
grid on;


ItFigOut = ItFigIn + 1;
fprintf('CPU time --> %.2f [min]\n\n',(cputime-TStart)/60);
end