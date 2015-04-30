function Pf = peng2004_wedge(X,P,v1,v2,steps,h)
% PENG2004_WEDGE moves a point with initial coordinates
% X towards the wedge with corner at P, and defined by vectors
% v1 and v2, p=3, as described in
% "Thick Surfaces - PhD Thesis", Peng [2004].
%
% Pf = peng2004_wedge(X,P,v1,v2,steps,h)
%
% Inputs:
%   X  #moving_points by 3 ccordinates of the initial moving points
%   P  3x1 wedge corner
%   v1  3x1 first vector of the wedge
%   v2  3x1 second vector of the wedge
%   steps   number of steps to be taken in the gradient flow
%   h    time step
% Outputs:
%   Pf   #moving_points by 3 matrix with final positions

for k=1:steps
    
    [int,grad_int] = peng2004_wedge_integral_gradient(X,P,v1,v2);
    
    grad = -[int.^(-1/3-1) int.^(-1/3-1) int.^(-1/3-1)].*grad_int;
    grad = normalizerow(grad);
    
    X = X+h*(-grad);
    input('');
    
end

Pf = X;