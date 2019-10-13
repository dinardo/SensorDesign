%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute particle's spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WorkTransportTotal = Total "Work-Transport" matrix
% x, y       = Axes
% NParticles = Total number of particles to be simulated
% Pitch      = Strip pitch [um]
% Bulk       = Bulk thickness [um]
% Radius     = Unit step of the movements and field interpolation [um]
% PType      = Particle type ['alpha' 'beta' 'gamma']
% ItFigIn    = Figure iterator input

function [ItFigOut] = ComputeSpectra(WorkTransportTotal,x,y,...
    NParticles,Pitch,Bulk,Radius,PType,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
eNoise  = 1000;  % Electronic noise [electrons]
nBins   = 100;   % Spectrum's number of bins
IsGamma = false; % Define whether the particle is a gamma or not

if strcmp(PType,'beta') == true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Charge spectrum for Beta %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ehLength = 70;            % Electron-holes pairs per unit length [electrons/um]
    depth    = Bulk;          % Source penetration depth [um]
    eMax     = ehLength*Bulk; % Maximum released charge [electrons]
elseif strcmp(PType,'alpha') == true || strcmp(PType,'gamma') == true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Charge spectrum for Alpha/Gamma %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    energyPair = 3.2;                  % Energy to create an electron-hole
                                       % pair [eV] [3.2 Silicon, 13 Diamond]
    eMax       = 4520000 / energyPair; % Maximum released charge [electrons]
    mean       = eMax;                 % Particle spectrum mean [electrons]
    sigma      = 17000 / energyPair;   % Particle spectrum sigma [electrons]
    % Particle penetration depth [um]
    if strcmp(PType,'alpha') == true
        depth = 10;
    elseif strcmp(PType,'gamma') == true
        depth   = Bulk * rand(1,1);
        IsGamma = true;
    end
else
    fprintf('Unknown particle: %s\n',PType);
    return;
end

eMax = eMax * 2;
EnergyScale = 0:eMax/nBins:eMax; % Spectrum energy axis [electrons]


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m calculating the spectrum of %d %s particles @@@\n',NParticles,PType);
HistoCharge = zeros(NParticles,1);

for i = 1:NParticles

    if strcmp(PType,'beta') == true
    % Particles enter randomly between -Pitch/2 -- +Pitch/2
    % Particles exit  randomly between -Pitch/2 -- +Pitch/2
%        enter = Pitch/2 * (2*rand(1,1) - 1); % x-coordinate entering particle
%        exit  = Pitch/2 * (2*rand(1,1) - 1); % x-coordinate exiting  particle
% @TMP@
        enter = 0;
        exit  = enter;
        mean  = ehLength*sqrt(depth^2 + (exit-enter)^2); % Landau MPV [electrons]
        sigma = mean / 12; % Scale factor between MPV and sigma of Landau [electrons]
        ChargeDensity = LandauRND(mean,sigma);
    elseif strcmp(PType,'alpha') == true || strcmp(PType,'gamma') == true
    % Particles enter randomly between -Pitch/2 -- +Pitch/2
    % Particles exit  randomly at Xin + Depth*tan(theta)
    % where theta is chosen randomly between -pi/4 and pi/4
        enter = 0;
        exit  = enter;
%        enter = Pitch/2 * (2*rand(1,1) - 1);
%        exit  = enter + depth*tan(pi/2*rand(1,1) - pi/4);
        ChargeDensity = normrnd(mean,sigma);
    end
    
    Noise         = normrnd(0,eNoise);
    ChargeDensity = ChargeDensity + Noise;
    if strcmp(PType,'gamma') == false
        ChargeDensity = ChargeDensity / sqrt(depth^2 + (exit-enter)^2);
    else
        ChargeDensity = ChargeDensity / Radius;
    end
    
    [Charge] = ComputeSignal(WorkTransportTotal,x,y,depth,...
        enter,exit,Bulk,Radius,ChargeDensity,IsGamma);
    
    HistoCharge(i) = Charge;
end


%%%%%%%%%
% Plots %
%%%%%%%%%
figure(ItFigIn);
histo = hist(HistoCharge,EnergyScale);
plot(EnergyScale,histo,'-o','LineWidth',1,...
    'MarkerEdgeColor','black','MarkerFaceColor','blue','MarkerSize',7);


%%%%%%%
% Fit %
%%%%%%%
if strcmp(PType,'beta') == true
    ft = fittype('LandauHist(x,mpv,sigma,ampl)');
    f = fit(EnergyScale',histo',ft,'StartPoint',[mean sigma NParticles/10]);
elseif strcmp(PType,'alpha') == true || strcmp(PType,'gamma') == true
    f = fit(EnergyScale',histo','gauss1');
end
hold on;
plot(f);
legend('Simulation','Fit');
hold off;

title('Collected charge histogram');
xlabel('Energy [e-]');
ylabel('Entries [#]');
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print fit results on screen %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Fit parameters:');
f


%%%%%%%%%%%%%%%%%%%%%%%
% Save data into file %
%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen('Spectrum.txt','w');
fprintf(fileID,'%12s %12s %12s %12s\n','Energy [e-]','Entries [#]','x-error','y-error');
xErrors = eMax/nBins*ones(size(EnergyScale));
yErrors = sqrt(histo);
for i = 1:length(histo)
    if yErrors(i) == 0
        xErrors(i) = 0;
    end
end
fprintf(fileID,'%12.2f %12.2f %12.2f %12.2f\n',[EnergyScale; histo;...
    xErrors; yErrors]);
fclose(fileID);

ItFigOut = ItFigIn + 1;
fprintf('CPU time --> %.2f [min]\n\n',(cputime-TStart)/60);
end