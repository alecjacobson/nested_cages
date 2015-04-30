function [int,grad_int,dist] = peng2004_canonical_wedge_integral_gradient(u,v,w,beta)

beta = beta*ones(size(u));
q = normrow([u v w]);

% if (abs(w)>1e-16)
    int = (2./w).*(atan2(w,(q-u).*cot(beta/2)-v));
% else
%     int = 2./((q-u).*cot(beta/2)-v);
% end
dist = int.^(-1/3);
% calculate gardient of the integral
grad_int(:,1) = sin(beta)./(q.*(q-(u.*cos(beta)+v.*sin(beta))));
grad_int(:,2) = 1./(q.*(q-u)) - cos(beta)./(q.*(q-(u.*cos(beta)+v.*sin(beta))));
% if (abs(w)>1e-16)
    grad_int(:,3) = -int./w +...
        ((u.^2+v.^2 - q.*u).*sin(beta)-q.*v.*(1-cos(beta)))./...
        (q.*w.*(q-u).*(q-(u.*cos(beta)+v.*sin(beta))));
% else
%     grad_int(:,3) = zeros(size(w));
% end

end