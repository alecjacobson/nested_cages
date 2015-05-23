function Pf = peng2004_surface(X,V0,F0,steps,h)
% PENG2004_SURFACE moves a point with initial coordinates
% X towards closed mesh, p=3, as described in
% "Thick Surfaces - PhD Thesis", Peng [2004].
%
% Pf = peng2004_surface(X,V0,F0,steps,h)
%
% Inputs:
%   X  #moving_points by 3 ccordinates of the initial moving points
%   (V0,F0)  mesh that defines the distance field
%   steps   number of steps to be taken in the gradient flow
%   h    time step
% Outputs:
%   Pf   #moving_points by 3 matrix with final positions

cameratoolbar;
trisurf(F0,V0(:,1),V0(:,2),V0(:,3),'FaceAlpha',0.1,'FaceColor',[0.5 0.0 0.0])
axis equal

for k=1:steps
    
    k
    
    % beggining of the C++ code (20x speedup)
    tic
    [int,grad_int] = peng2004_surface_integral_gradient_mex(X,V0,F0)
    toc
    
%     tic
%     
%     int = zeros(size(X,1),1);
%     grad_int = zeros(size(X,1),3);
%     
%     for j=1:size(F0,1)
%     
%         v1 = V0(F0(j,1),:)';
%         v2 = V0(F0(j,2),:)';
%         v3 = V0(F0(j,3),:)';
%         % First wedge
%         P_1 = v1; e1_1 = v2-v1; e2_1 = v3-v1;
%         % Second wedge
%         P_2 = v3; e1_2 = v2-v3; e2_2 = v3-v1;
%         % Third wedge
%         P_3 = v2; e1_3 = v2-v3; e2_3 = v2-v1;
% 
%         [int_1,grad_int_1] = peng2004_wedge_integral_gradient(X,P_1,e1_1,e2_1);
%         [int_2,grad_int_2] = peng2004_wedge_integral_gradient(X,P_2,e1_2,e2_2);
%         [int_3,grad_int_3] = peng2004_wedge_integral_gradient(X,P_3,e1_3,e2_3);
% 
%         int = int + (int_1-int_2+int_3);
%         grad_int = grad_int + (grad_int_1-grad_int_2+grad_int_3);
%         
%     end
%     
%     int
%     grad_int
%     
%     toc
    
    % end of the C++ code
    
    grad = -[int.^(-1/3-1) int.^(-1/3-1) int.^(-1/3-1)].*grad_int;
    grad = normalizerow(grad);
    signs = -2*(winding_number(V0,F0,X))+1;
    grad = grad.*[signs signs signs];
    
    X = X+h*(-grad)
    hold on
    plot3(X(:,1),X(:,2),X(:,3),'.k','MarkerSize',10);
    drawnow
    input('');
    
end

Pf = X;