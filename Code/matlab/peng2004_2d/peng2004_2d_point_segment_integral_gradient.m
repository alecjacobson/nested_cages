function [grad,int,grad_int,distance] = peng2004_2d_point_segment_integral_gradient(X,a0,b0,a1,b1,a)
% Calculates the integral that defines the p-distance of a point (x,y) 
% to a segment (a0,b0)-(a1,b1). Rotates the coordinate system and 
% calculates the gradient and the integral w.r.t. [0,A]
%
% [grad,int,grad_int,distance] = peng2004_2d_point_segment_integral_gradient(X,a0,b0,a1,b1,a)
%
% Input:
%   X:  #moving_points by 1 x-ccordinates of the initial moving points
%   (a0,b0)-(a1,b1):  endpoints of the segment
%   a:    exponent that defines the p-distance in Peng et al. [2004]
% Output:
%   grad:   #moving_points by 2 full gradient
%   int:  #moving_points by 1 integral value
%   grad:   #moving_points by 2 full gradient of the integral (only)
%   distance: #moving_points by 1 full distance

if (a1==a0 && b1 > b0)
    theta = -pi/2;
elseif (a1==a0 && b1 <= b0)
    theta = pi/2;
else
    theta = -atan2((b1-b0),(a1-a0));
end
% Make a function grad_segment for this repeated things
A = sqrt((a1-a0)^2+(b1-b0)^2)*ones(size(X,1),1);
a0 = a0*ones(size(X,1),1);
b0 = b0*ones(size(X,1),1);
X_tilde = ([cos(theta) -sin(theta); sin(theta) cos(theta)]*[X(:,1)-a0 X(:,2)-b0]')';
int = peng2004_2d_point_segment_0A_integral(X_tilde(:,1),X_tilde(:,2),A,a);
grad_int = peng2004_2d_point_segment_0A_integral_gradient(X_tilde(:,1),X_tilde(:,2),A,a);
distance = int.^(-1/a);
grad = -[int.^(-1/a-1) int.^(-1/a-1)].*grad_int;
grad = normalizerow(grad);
grad = ([cos(-theta) -sin(-theta); sin(-theta) cos(-theta)]*grad')';
% have to rotate grad_int to output it
grad_int = ([cos(-theta) -sin(-theta); sin(-theta) cos(-theta)]*grad_int')';

end