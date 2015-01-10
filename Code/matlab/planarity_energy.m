  function [E,cb_data] = planarity_energy(V,F_quad)
    % Input:
    %  V: input triangle mesh
    %  F_quad: #quadsx4 list of vertex indices in V. 
    %     The energy will measure the square sum of coplanarity of each set of 
    %     vertices V(F_quad(i,1),:),V(F_quad(i,2),:),V(F_quad(i,3),:),V(F_quad(i,4),:)
    %  Output:
    %  E: energy value
    
    vol = zeros(size(F_quad,1),1);
    for k=1:size(F_quad)
        vol(k) = det([V(F_quad(k,1),:) 1;V(F_quad(k,2),:) 1; ...
            V(F_quad(k,3),:) 1; V(F_quad(k,4),:) 1]);
    end
    E = sum(vol.^2);
    cb_data = [];
    
  end