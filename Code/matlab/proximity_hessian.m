function H = proximity_hessian(P,P_coarse,u,inds)
% PROXIMITY_HESSIAN
% H = proximity_hessian(P,P_coarse,u,inds)
%
% Calculate the hessian to minimize proximity energy
%
% Input:
%   P  2xn_all list of polygon vertex positions of the fine mesh
%   P_coarse 2xn_coarse list of polygon vertex positions of the fine mesh
%   u:  1xn_all list orthogonal projection coefficients for each point in P
%   inds: indices of the points in P_coarse in P's order
% Output:
%   H: (n_all+n+coarse)x(n_all+n+coarse) hessian matrix

P = P';
P_coarse = P_coarse';
H = sparse(size(P_coarse,1),size(P_coarse,1));

for k=1:size(inds,2)-1
        
    for j=inds(k):inds(k+1)-1
        
        H(k,k) = H(k,k)+2*(1-u(j))^2;
        H(k,k+1) = H(k,k+1)+2*(1-u(j))*u(j);
        H(k+1,k+1) = H(k+1,k+1)+2*u(j)^2;
        
    end
    
end

% force is negative of the gradient
H = blkdiag(-H,sparse(size(P,1),size(P,1)));