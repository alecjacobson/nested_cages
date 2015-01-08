  function [grad,cb_data] = planarity_gradient(V,F_quad)
    % Input:
    %  V: input triangle mesh
    %  F_quad: #quadsx4 list of vertex indices in V. 
    %     The energy will measure the square sum of coplanarity of each set of 
    %     vertices V(F_quad(i,1),:),V(F_quad(i,2),:),V(F_quad(i,3),:),V(F_quad(i,4),:)
    %  Output:
    %  grad: energy gradient
    
    grad = zeros(size(V));
    minors = zeros(4*size(F_quad,1),1);
    for k=1:size(F_quad)
        A = [V(F_quad(k,1),:) 1;V(F_quad(k,2),:) 1; ...
            V(F_quad(k,3),:) 1; V(F_quad(k,4),:) 1];
        vol = det(A);
%         for i=1:4
%             temp = A;
%             for j=1:3
%                 temp(i,:) = [];
%                 temp(:,j) = [];
%                 minors(4*(k-1)+i,j) = (-1)^(i+j)*det(temp);
%                 temp = A;
%             end
%         end
        temp = A;
        for i=1:4
            for j=1:3
                v = F_quad(k,i);
                temp(i,:) = [];
                temp(:,j) = [];
                grad(v,j) = grad(v,j) + vol*((-1)^(i+j))*det(temp);
                temp = A;
            end
        end
               
        
    end
    
    cb_data = [];
    
  end
  
% Test that works for step+project  
%   for k=1:2000
% grad = planarity_gradient(V,F_test);
% [V_eltopo,rest_dt] = collide_eltopo_mex([V0;V],[F0;size(V0,1)+cages_F{1}],[V0;V-(5e-1)*grad],size(V0,1),1e-3,1e-1);
% V = V_eltopo(size(V0,1)+1:end,:);
% energy = planarity_energy(V,F_test)
% k
% end