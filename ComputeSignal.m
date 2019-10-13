%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the signal given the "Work-Transport" matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WorkTransportTotal = Total "Work-Transport" matrix
% x, y           = Axes
% Depth          = Penetration depth along y (positive = from back side,
%                  negative = from strip side) [um]
% XEnterParticle = X cordinate of the enter of the particle [um] (back-plane side)
% XExitParticle  = X cordinate of the exit of the particle [um] (strip side)
% Bulk           = Bulk thickness [um]
% Radius         = Unit step of the movements and field interpolation [um]
% ChargeDensity  = Charge density [electrons/um]
% IsGamma        = [true/false]

function [Charge] = ComputeSignal(WorkTransportTotal,x,y,Depth,...
    XEnterParticle,XExitParticle,Bulk,Radius,ChargeDensity,IsGamma)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Charge = 0;         % Starting value of the signal [electrons]
eps    = Radius/10; % Minimum step to avoid infinities [um]


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
mx = min(abs(x - XEnterParticle));
if Depth >= 0
    my = y(1);
else
    my = min(abs(y - (Bulk + Depth)));
end

while my >= 0 && ((Depth >= 0 && my <= Depth) || (Depth < 0 && my <= Bulk))

    x_near = find(abs(x - mx) <= Radius);
    y_near = find(abs(y - my) <= Radius);

    % Calculate the average Work within a radius = Radius
    W      = 0;
    weight = 0;
    for ii = 1:length(x_near)
        for jj = 1:length(y_near)

            dist = sqrt((mx-x(x_near(ii)))^2 + (my-y(y_near(jj)))^2);

            if dist <= Radius && ~isnan(WorkTransportTotal(y_near(jj),x_near(ii)))
                
                if dist < eps
                    dist = eps;
                end
                
                W = W + WorkTransportTotal(y_near(jj),x_near(ii)) / dist;
                weight = weight + 1 / dist;
            end
        end
    end
    if weight ~= 0
        W = W / weight;
    end
    
    Charge = Charge + W;

    % Calculate the next movement
    if IsGamma == false
        dx   = XExitParticle - XEnterParticle;
        sinv = dx   / sqrt(dx^2 + Bulk^2);
        cosv = Bulk / sqrt(dx^2 + Bulk^2);
        mx   = mx + Radius*sinv;
        my   = my + Radius*cosv;
    else
        break
    end
end

Charge = ChargeDensity * Charge;
end
