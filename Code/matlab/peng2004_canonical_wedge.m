function Pf = peng2004_canonical_wedge(X,beta,steps,h)
% PENG2004_CANONICAL_WEDGE moves a point with initial coordinates
% X towards the wedge with corner at (0,0,0), one of the vectors
% pointing towards the positive x direction and the other making angle
% beta with the x-axis (0 <= beta <= pi), p=3, as described in
% "Thick Surfaces - PhD Thesis", Peng [2004].
%
% Pf = peng2004_canonical_wedge(X,beta,steps,h)
%
% Inputs:
%   X  #moving_points by 3 ccordinates of the initial moving points
%   beta  wedge angle
%   steps   number of steps to be taken in the gradient flow
%   h    time step
% Outputs:
%   Pf   #moving_points by 3 matrix with final positions

u = X(:,1);
v = X(:,2);
w = X(:,3);

for k=1:steps
    
    [int,grad_int,dist] = peng2004_canonical_wedge_integral_gradient(u,v,w,beta);
    
    int
    grad_int
    dist
    
    grad(:,1) = -(dist.^4).*grad_int(:,1);
    grad(:,2) = -(dist.^4).*grad_int(:,2);
    grad(:,3) = -(dist.^4).*grad_int(:,3);

%     grad = -[int.^(-1/3-1) int.^(-1/3-1) int.^(-1/3-1)].*grad_int;
%     grad
    
    grad = normalizerow(grad)
    
    X = real(X+h*(-grad))
    u = X(:,1);
    v = X(:,2);
    w = X(:,3);
%     input('');
    
    
end

Pf = X;