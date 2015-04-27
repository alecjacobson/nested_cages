function grad_I = peng2004_2d_point_segment_0A_integral_gradient(x,y,A,a)
% Calculates the gradient of the integral that defines the p-distance of
% a point (x,y) to a segment [0,A] (y=0). Details in
% "Interactively Modeling of Topologically Complex Geometric Detail"
% Peng et a. [2004].
%
% grad_I = peng2004_2d_point_segment_0A_integral_gradient(x,y,A,a)
%
% Input:
%   x:  #moving_points by 1 x-ccordinates of the initial moving points
%   y:  #moving_points by 1 y-ccordinates of the initial moving points
%   A:  surface is the segment [0,A]
%   a:    exponent that defines the p-distance in Peng et al. [2004]
% Output:
%   Pf:   #moving_points by 2 vector with final positions
        
        if (a==1)
            grad_I(:,1) = (x.^2+y.^2).^(-1/2)+(-1).*(A.^2+(-2).*A.*x+x.^2+y.^2).^(-1/2);
            grad_I(:,2) = ((-1).*A+x+((A+(-1).*x).^2+y.^2).^(1/2)).*(x+(x.^2+y.^2).^(1/2)) ...
                .^(-1).*(y.*(x.^2+y.^2).^(-1/2).*((-1).*A+x+((A+(-1).*x).^2+y.^2) ...
                .^(1/2)).^(-1)+(-1).*y.*((A+(-1).*x).^2+y.^2).^(-1/2).*((-1).*A+x+ ...
                ((A+(-1).*x).^2+y.^2).^(1/2)).^(-2).*(x+(x.^2+y.^2).^(1/2)));
        elseif (a==2)
            grad_I(:,1) = A.*(A+(-2).*x).*(x.^2+y.^2).^(-1).*(A.^2+(-2).*A.*x+x.^2+y.^2).^( ...
                -1);
            grad_I(:,2) = (-1).*A.*y.^(-1).*(A.*x+(-1).*x.^2+y.^2).*(x.^2+y.^2).^(-1).*( ...
                A.^2+(-2).*A.*x+x.^2+y.^2).^(-1)+(-1).*y.^(-2).*(atan((A+(-1).*x)./ ...
                y)+atan(x./y));
        elseif (a==3)
            grad_I(:,1) = (x.^2+y.^2).^(-3/2)+(-1).*(A+(-1).*x).*y.^(-2).*((A+(-1).*x).^2+ ...
                y.^2).^(-3/2).*((-1).*A+x+((A+(-1).*x).^2+y.^2).^(1/2))+(-1).*y.^( ...
                -2).*(A.^2+(-2).*A.*x+x.^2+y.^2).^(-1).*((-1).*A+x+(A.^2+(-2).*A.* ...
                x+x.^2+y.^2).^(1/2));
            grad_I(:,2) = y.^(-3).*((-1).*A.*(A.^2+(-2).*A.*x+x.^2+y.^2).^(-3/2).*(2.*A.^2+( ...
                -4).*A.*x+2.*x.^2+3.*y.^2)+x.*((-2).*(x.^2+y.^2).^(-1/2)+2.*(A.^2+ ...
                (-2).*A.*x+x.^2+y.^2).^(-1/2)+y.^2.*((-1).*(x.^2+y.^2).^(-3/2)+( ...
                A.^2+(-2).*A.*x+x.^2+y.^2).^(-3/2))));
        else
            error('Implemented only for a=1,2 or 3')
        end
            
    end