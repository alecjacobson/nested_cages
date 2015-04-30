function Pf = peng2004_triangle(X,v1,v2,v3,steps,h)
% PENG2004_TRIANGLE moves a point with initial coordinates
% X towards a tirnagle with vertices v1,v2,v3, p=3, as described in
% "Thick Surfaces - PhD Thesis", Peng [2004].
%
% Pf = peng2004_triangle(X,v1,v2,v3,steps,h)
%
% Inputs:
%   X  #moving_points by 3 ccordinates of the initial moving points
%   P  3x1 wedge corner
%   v1  3x1 first vertex of the triangle
%   v2  3x1 second vertex of the triangle
%   v3  3x1 third vertex of the triangle
%   steps   number of steps to be taken in the gradient flow
%   h    time step
% Outputs:
%   Pf   #moving_points by 3 matrix with final positions

% First wedge
P_1 = v1; e1_1 = v2-v1; e2_1 = v3-v1;
% Second wedge
P_2 = v3; e1_2 = v2-v3; e2_2 = v3-v1;
% Third wedge
P_3 = v2; e1_3 = v2-v3; e2_3 = v2-v1;
cameratoolbar;
trisurf([1 2 3],[v1(1);v2(1);v3(1)],[v1(2);v2(2);v3(2)],[v1(3);v2(3);v3(3)],'FaceAlpha',0.3)

for k=1:steps
    
    [int_1,grad_int_1] = peng2004_wedge_integral_gradient(X,P_1,e1_1,e2_1);
    [int_2,grad_int_2] = peng2004_wedge_integral_gradient(X,P_2,e1_2,e2_2);
    [int_3,grad_int_3] = peng2004_wedge_integral_gradient(X,P_3,e1_3,e2_3);
    
    int = int_1-int_2+int_3;
    grad_int = grad_int_1-grad_int_2+grad_int_3;
    
    grad = -[int.^(-1/3-1) int.^(-1/3-1) int.^(-1/3-1)].*grad_int;
    grad = normalizerow(grad);
    
    X = X+h*(-grad)
    hold on
    plot3(X(:,1),X(:,2),X(:,3),'.k','MarkerSize',10);
    input('');
    
end

Pf = X;