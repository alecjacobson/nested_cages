function I = peng2004_2d_point_segment_0A_integral(x,y,A,a)
% Calculates the integral that defines the p-distance of a point (x,y) 
% to a segment [0,A] (y=0). Details in
% "Interactively Modeling of Topologically Complex Geometric Detail"
% Peng et a. [2004].
%
% I = peng2004_2d_point_segment_0A_integral(x,y,A,a)
%
% Input:
%   x:  #moving_points by 1 x-ccordinates of the initial moving points
%   y:  #moving_points by 1 y-ccordinates of the initial moving points
%   A:  surface is the segment [0,A]
%   a:    exponent that defines the p-distance in Peng et al. [2004]
% Output:
%   Pf:   #moving_points by 2 vector with final positions

if (a==1)
    I = log(((-1).*A+x+((A+(-1).*x).^2+y.^2).^(1/2)).^(-1).*(x+(x.^2+y.^2) ...
        .^(1/2)));
elseif (a==2)
    I = y.^(-1).*(atan((A+(-1).*x)./y)+atan(x./y));
    % changed to atan2
elseif (a==3)
    I = 2.*y.^(-2)+(-1).*y.^(-2).*((A+(-1).*x).^2+y.^2).^(-1/2).*((-1).*A+ ...
        x+((A+(-1).*x).^2+y.^2).^(1/2))+(-1).*y.^(-2).*(1+(-1).*x.*(x.^2+ ...
        y.^2).^(-1/2));
else
    error('Implemented only for a=1,2 or 3')
end

end