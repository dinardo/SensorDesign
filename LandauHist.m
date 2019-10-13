%%%%%%%%%%%%%%%%%%%%
% Landau histogram %
%%%%%%%%%%%%%%%%%%%%

function y = LandauHist(x,mpv,sigma,ampl)

y = zeros(size(x));
for i = 1:length(x)
    y(i) = ampl * Landau(x(i),mpv,sigma,true);
end
end
