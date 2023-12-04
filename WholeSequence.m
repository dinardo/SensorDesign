%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO                                                     %
% - Review transport equations in magnetic field            %
% - Define Sensor&Air volumes in SolvePoisson3D_PlanarPixel %
%   (not available in MATLAB 2016)                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Clean up everything
close all;
clear;
clc;


         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Program to caculate the signal in particle detectors %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
BiasV = -120; % Sensor backplane voltage [V]

Fluence = 0; % Irradiation fluence [10^16 1MeV n.eq./cm^2]
             % 1/tau = c*Fluence/(1 + c*Fluence/t), extracted from fit to data [ns^-1]
ce = 5.36;
te = 0.8295;
ch = 3.361;
th = 107.6;
scale = 2.7; % Scale factor to correct the data life-time [tuned on ATLAS results]
TauBe = scale * (1 + ce*Fluence/te)/(ce*Fluence); % Life-time on the backplane side [ns]
TauSe = scale * (1 + ce*Fluence/te)/(ce*Fluence); % Life-time on the strip side [ns]
TauBh = scale * (1 + ch*Fluence/th)/(ch*Fluence); % Life-time on the backplane side [ns]
TauSh = scale * (1 + ch*Fluence/th)/(ch*Fluence); % Life-time on the strip side [ns]

Bulk   =   200; % Bulk thickness [um]
PitchX =   250; % Pitch along X [um]
PitchY = 20000; % Pitch along Y [um]

qe       = -1.6e-19; % Electron charge [Coulomb]
eps0     = 8.85e-18; % Vacuum permittivity [F/um]
epsR     = 11.7;     % Relative permittivity [11.7 Silicon, 5.7 Diamond, 12.85 GaAs]
dN_dPhi  = 35;       % dN/dPhi extracted from data [#/(um^3 10^16)]
DeplVnoF = 10;       % Full depletion voltage for non irradiated sensors [V]
DeplV    = qe*Bulk^2/(2*epsR*eps0)*dN_dPhi*Fluence - DeplVnoF; % Sensor full depletion voltage [V]
rho      = 2*DeplV*epsR*eps0/(qe*Bulk^2); % Bulk doping concentration Planar Strip [#/um^3]
%rho      = 2*DeplV*epsR*eps0/(qe*sqrt(PitchX^2+PitchY^2)/2); % Bulk doping concentration 3D Pixel [#/um^3]

BField = 0.0; % Magnetic field (orthogonal+outgoing from 2D geometry) [T]

T = 300; % Sensor temperature [Kelvin]

mu_e   = 140*(T/300)^(-2.4); % Electron mobility [um^2/(V*ns)] [140 Silicon, 180 Diamond, 850 GaAs]
RH_e   = 1;   % Relative Hall electron mobility [1 Silicon, 1 Diamond]
vs_e   = 110; % Saturation velocity of the electrons [um/ns] [110 Silicon, 260 Diamond, 200 GaAs]
beta_e = 0.0257*T^0.66; % Exponent for the electric field dependence of the mobility

mu_h   = 48*(T/300)^(-2.2); % Hole mobility in [um^2/(V*ns)] [45 Silicon, 120 Diamond, 45 AsGa]
RH_h   = 1;    % Relative Hall hole mobility in [1 Silicon, 1 Diamond] 
vs_h   = 95;   % Saturation velocity of the holes [um/ns] [95 Silicon, 160 Diamond, 95 GaAs]
beta_h = 0.46*T^0.17; % Exponent for the electric field dependence of the mobility

Step   = 2;       % Unit step of the lattice on which the field is computed [um]
Radius = Step/10; % Unit step of the movements and field interpolation [um]

XQ = 0; % Coordinate for potential query along z [um]
YQ = 0; % Coordinate for potential query along z [um]

NAverage   = 2;      % Generate NAverage "Work-Transport" matrices and average them
NParticles = 100;    % Total number of particles to be simulated
PType      = 'beta'; % Particle type ['alpha' 'beta' 'gamma']

fprintf('@@@ Derived parameters @@@\n');
fprintf('\t- Electron''s mobility --> %.1f(T=%.1f) [um^2/(V ns)]\n',mu_e,T);
fprintf('\t- Hole''s mobility --> %.1f(T=%.1f) [um^2/(V ns)]\n',mu_h,T);
fprintf('\t- Electron''s life-time --> %.2f [ns], %.2f [ns]\n',TauBe,TauSe);
fprintf('\t- Hole''s life-time --> %.2f [ns], %.2f [ns]\n',TauBh,TauSh);
fprintf('\t- Full depletion voltage --> %.1f [V]\n',DeplV);
fprintf('\t- Depleted depth --> %.1f [um]\n',Bulk*sqrt(BiasV/DeplV));
fprintf('\t- Doping concentration --> %.1E [#/cm^3]\n',rho*1e12);
fprintf('\t- Resistivity --> %.1E [Ohm cm]\n\n',-1/(qe*mu_h*rho)*1e-13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


do2D3Dcheck = false;
doSignal    = false;
rng default; % Reset random seed
ItFig = 1;   % Figure iterator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recompute sensor thickness based on the depth of depleted region %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Bulk * sqrt(BiasV/DeplV) < Bulk
    Bulk = Bulk * sqrt(BiasV/DeplV);
    fprintf('@@@ Changed bulk thickness to %.1f [um] @@@\n\n',Bulk);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the potentials %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Planar Strip
[pdem, TotalPot, DecomposedGeom, BulkStart, BulkStop, VolumeHeight] = StripPlanar_SolvePoisson2D(Bulk,PitchX,BiasV,epsR,rho*qe);
[~, ~, ItFig] = Planar_Plots(pdem,TotalPot,DecomposedGeom,BulkStart,BulkStop,VolumeHeight,PitchX,PitchY,0,epsR,true,XQ,ItFig);

[pdem, WeightPot, DecomposedGeom, BulkStart, BulkStop, VolumeHeight] = StripPlanar_SolvePoisson2D(Bulk,PitchX,0,epsR,0);
[Sq2D, xq2D, ItFig] = Planar_Plots(pdem,WeightPot,DecomposedGeom,BulkStart,BulkStop,VolumeHeight,PitchX,PitchY,1,epsR,true,XQ,ItFig);

% 3D Pixel
% [pdem, TotalPot, DecomposedGeom] = Pixel3D_SolvePoisson2D(PitchX,PitchY,BiasV,epsR,rho*qe);
% [~, ~, ItFig] = Planar_Plots(pdem,TotalPot,DecomposedGeom,0,Bulk,0,PitchX,PitchY,0,epsR,false,XQ,ItFig);
% 
% [pdem, WeightPot, DecomposedGeom] = Pixel3D_SolvePoisson2D(PitchX,PitchY,0,epsR,0);
% [~, ~, ItFig] = Planar_Plots(pdem,WeightPot,DecomposedGeom,0,Bulk,0,PitchX,PitchY,1,epsR,false,XQ,ItFig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the potential in 2D and 3D %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do2D3Dcheck == true
    fprintf('@@@ I''m comparing the potential for 2D and 3D @@@\n');
    [~, Sq3D, xq3D, ItFig] = PlanarPixel_SolvePoisson3D(Bulk,PitchX,PitchY,0,1,epsR,0,XQ,YQ,ItFig);
    Diff2D3D = ((Sq2D ./ xq2D' + Sq2D ./ (Bulk - xq2D')) ./...
        (Sq3D ./ xq3D' + Sq3D ./ (Bulk - xq3D')) - 1) * 100;
    fprintf('@@@ Potential percentage difference %.1f%% @@@\n\n',mean(Diff2D3D,'omitnan'));

    figure(ItFig);
    plot(xq2D,Diff2D3D);
    title(sprintf('Potential percentage difference at x = %.2f um y = %.2f um',XQ,YQ));
    xlabel('Z [\mum]');
    ylabel('Percentage [%]');
    grid on;
    ItFig = ItFig + 1;
end


if doSignal == true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the velocity field %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [VFieldx_e, VFieldy_e, VFieldx_h, VFieldy_h, x, y, ItFig] =...
        VelocityField(TotalPot,Step,Bulk,BField,PitchX,...
        mu_e,RH_e,vs_e,beta_e,mu_h,RH_h,vs_h,beta_h,ItFig);


    %%%%%%%%%%%%%%%%%%%%
    % Compute the work %
    %%%%%%%%%%%%%%%%%%%%
    [WorkTransportTotal, x, y, ItFig] =...
        ManyWorkTransport(WeightPot,VFieldx_e,VFieldy_e,VFieldx_h,VFieldy_h,...
        x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh,NAverage,ItFig);


    %%%%%%%%%%%%%%%%%
    % Fit/Smoothing %
    %%%%%%%%%%%%%%%%%
    fprintf('@@@ I''m interpolating the work transport matrix with a step of %.2f um @@@\n\n',Radius);
    subx = x(1):Radius:x(length(x));
    suby = y(1):Radius:y(length(y));
    [ssubx,ssuby] = meshgrid(subx,suby);
    [x0, y0, z0] = prepareSurfaceData(x,y,WorkTransportTotal);
    [WorkTransportTotal_, goodness, output] = fit([x0 y0],z0,'cubicinterp');
    subWorkTransportTotal = WorkTransportTotal_(ssubx, ssuby);

    figure(ItFig);
    colormap jet;
    surf(ssubx,ssuby,subWorkTransportTotal,'EdgeColor','none','FaceColor','interp');
    title('Interpolated Total <Work-Transport>');
    xlabel('X [\mum]');
    ylabel('Z [\mum]');
    zlabel('Work / q [#charges * V]');
    ItFig = ItFig + 1;


    %%%%%%%%%%%%%%%%%%%%%%
    % Compute the spectra %
    %%%%%%%%%%%%%%%%%%%%%%%
    [ItFig] = ComputeSpectra(subWorkTransportTotal,subx,suby,NParticles,...
        PitchX,Bulk,Radius,PType,ItFig);


    %%%%%%%%%%%%%%%%%%%%%%%%
    % Alert for completion %
    %%%%%%%%%%%%%%%%%%%%%%%%
    load handel;
    sound(y,Fs);
end
