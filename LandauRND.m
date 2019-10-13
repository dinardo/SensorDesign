%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Landau-distributed random numbers %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = LandauRND(mpv,sigma)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSigma = 40;
xAxis  = mpv + sigma*nSigma;
toss   = true;


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
while toss == true
    x    = rand(1,1)*xAxis;
    coin = rand(1,1)*0.180655; % Maximum of Landau

    if coin < Landau(x,mpv,sigma,false);
        out  = x;
        toss = false;
    end
end
end
