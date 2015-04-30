function [int,grad_int] = peng2004_wedge_integral_gradient(X,P,v1,v2)

% first rotation: to make v1(3) = 0;
theta_z = atan2(v1(3),v1(2));
R_theta_z = [1 0 0; 0 cos(-theta_z) -sin(-theta_z); 0 sin(-theta_z) cos(-theta_z)];
% translate and rotate moving points
X = (R_theta_z*[X(:,1)-P(1) X(:,2)-P(2) X(:,3)-P(3)]')';
% rotate the wedge vectors
v1 = R_theta_z*v1;
v2 = R_theta_z*v2;
% second rotation to make v1(2)=0
theta_y = atan2(v1(2),v1(1));
R_theta_y = [cos(-theta_y) -sin(-theta_y) 0; sin(-theta_y) cos(-theta_y) 0; 0 0 1];
X = (R_theta_y*[X(:,1) X(:,2) X(:,3)]')';
v1 = R_theta_y*v1;
v2 = R_theta_y*v2;
% third rotation to make v2(3) = 0
theta_3 = atan2(v2(3),v2(2));
R_theta_3 = [1 0 0; 0 cos(-theta_3) -sin(-theta_3); 0 sin(-theta_3) cos(-theta_3)];
X = (R_theta_3*[X(:,1) X(:,2) X(:,3)]')';
v1 = R_theta_3*v1;
v2 = R_theta_3*v2;
% define beta angle - to belong to [0,2*pi]
beta = atan2(v2(2),v2(1));
[int,grad_int,~] = peng2004_canonical_wedge_integral_gradient(X(:,1),X(:,2),X(:,3),beta);
% back to original coordinate system
grad_int = ([1 0 0; 0 cos(theta_3) -sin(theta_3); 0 sin(theta_3) cos(theta_3)]*grad_int')';
grad_int = ([cos(theta_y) -sin(theta_y) 0; sin(theta_y) cos(theta_y) 0; 0 0 1]*grad_int')';
grad_int = ([1 0 0; 0 cos(theta_z) -sin(theta_z); 0 sin(theta_z) cos(theta_z)]*grad_int')';