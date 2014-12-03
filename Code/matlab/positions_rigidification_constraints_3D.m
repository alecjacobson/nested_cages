function [A_rig,b_rig] = positions_rigidification_constraints_3D(EdgeCollisions,CV,EE_all)
  % POSITIONS_RIGIDIFICATION_CONSTRAINTS_3D
  %
  % positions_rigidification_constraints_3D(B,CV)
  %
  % Calculates matrix and rhs for ridifying primitives that are
  % intersecting and not being resolved on the coarse mesh
  %
  % Input:
  %   EdgeCollisions
  %   CV  positions of the vertices of the coarse mesh at previous
  %      (resolved) state
  %   EE_all edges of both fine and coarse mesh (w.r.t. EdgeCollisions is
  %   calculated)
  % Output:
  %   A_rig   (6*3*(number of collisions))x(3*num.vert. coarse mesh) matrices
  %           with entries +1 and -1
  %   b_rig   (6*3*(number of collisions))x1 rhs vector
  
  nun = size(CV,1);
  
  [EdgeColI,EdgeColJ,~] = find(EdgeCollisions);
  nc = max(max(EE_all))-nun;
  
  edge_num_collisions = size(EdgeColI,1);
  
  I_A_rig = ones(1,36*edge_num_collisions);
  J_A_rig = ones(1,36*edge_num_collisions);
  Val_A_rig = zeros(1,36*edge_num_collisions);
  
  b_rig = zeros(6*3*edge_num_collisions,1);
  
  for k=1:edge_num_collisions
      
      i = EdgeColI(k);
      j = EdgeColJ(k);
      
      if (i<=nc || j<=nc)
          error('one of the vertices belong to fine mesh. It should not be entering this function');
      end
      
      v1 = EE_all(i,1)-nc; w1 = EE_all(i,2)-nc;
      v2 = EE_all(j,1)-nc; w2 = EE_all(j,2)-nc;
      
      I_A_rig(36*(k-1)+1) = 6*3*(k-1)+1;
      J_A_rig(36*(k-1)+1) = v1;
      Val_A_rig(36*(k-1)+1) = 1;
      
      I_A_rig(36*(k-1)+2) = 6*3*(k-1)+1;
      J_A_rig(36*(k-1)+2) = w1;
      Val_A_rig(36*(k-1)+2) = -1;
      
      b_rig(6*3*(k-1)+1) = CV(v1,1) - CV(w1,1);
      
      
      I_A_rig(36*(k-1)+3) = 6*3*(k-1)+2;
      J_A_rig(36*(k-1)+3) = v1+nun;
      Val_A_rig(36*(k-1)+3) = 1;
      
      I_A_rig(36*(k-1)+4) = 6*3*(k-1)+2;
      J_A_rig(36*(k-1)+4) = w1+nun;
      Val_A_rig(36*(k-1)+4) = -1;
      
      b_rig(6*3*(k-1)+2) = CV(v1,2) - CV(w1,2);     
      
      I_A_rig(36*(k-1)+5) = 6*3*(k-1)+3;
      J_A_rig(36*(k-1)+5) = v1+2*nun;
      Val_A_rig(36*(k-1)+5) = 1;
      
      I_A_rig(36*(k-1)+6) = 6*3*(k-1)+3;
      J_A_rig(36*(k-1)+6) = w1+2*nun;
      Val_A_rig(36*(k-1)+6) = -1;
      
      b_rig(6*3*(k-1)+3) = CV(v1,3) - CV(w1,3);
      
      
      I_A_rig(36*(k-1)+7) = 6*3*(k-1)+4;
      J_A_rig(36*(k-1)+7) = v2;
      Val_A_rig(36*(k-1)+7) = 1;
      
      I_A_rig(36*(k-1)+8) = 6*3*(k-1)+4;
      J_A_rig(36*(k-1)+8) = w2;
      Val_A_rig(36*(k-1)+8) = -1;
      
      b_rig(6*3*(k-1)+4) = CV(v2,1) - CV(w2,1);
      
      I_A_rig(36*(k-1)+9) = 6*3*(k-1)+5;
      J_A_rig(36*(k-1)+9) = v2+nun;
      Val_A_rig(36*(k-1)+9) = 1;
      
      I_A_rig(36*(k-1)+10) = 6*3*(k-1)+5;
      J_A_rig(36*(k-1)+10) = w2+nun;
      Val_A_rig(36*(k-1)+10) = -1;
      
      b_rig(6*3*(k-1)+5) = CV(v2,2) - CV(w2,2);
      
      I_A_rig(36*(k-1)+11) = 6*3*(k-1)+6;
      J_A_rig(36*(k-1)+11) = v2+2*nun;
      Val_A_rig(36*(k-1)+11) = 1;
      
      I_A_rig(36*(k-1)+12) = 6*3*(k-1)+6;
      J_A_rig(36*(k-1)+12) = w2+2*nun;
      Val_A_rig(36*(k-1)+12) = -1;
      
      b_rig(6*3*(k-1)+6) = CV(v2,3) - CV(w2,3);

      
      I_A_rig(36*(k-1)+13) = 6*3*(k-1)+7;
      J_A_rig(36*(k-1)+13) = v1;
      Val_A_rig(36*(k-1)+13) = 1;
      
      I_A_rig(36*(k-1)+14) = 6*3*(k-1)+7;
      J_A_rig(36*(k-1)+14) = v2;
      Val_A_rig(36*(k-1)+14) = -1;
      
      b_rig(6*3*(k-1)+7) = CV(v1,1) - CV(v2,1);
      
      I_A_rig(36*(k-1)+15) = 6*3*(k-1)+8;
      J_A_rig(36*(k-1)+15) = v1+nun;
      Val_A_rig(36*(k-1)+15) = 1;
      
      I_A_rig(36*(k-1)+16) = 6*3*(k-1)+8;
      J_A_rig(36*(k-1)+16) = v2+nun;
      Val_A_rig(36*(k-1)+16) = -1;
      
      b_rig(6*3*(k-1)+8) = CV(v1,2) - CV(v2,2);
      
      I_A_rig(36*(k-1)+17) = 6*3*(k-1)+9;
      J_A_rig(36*(k-1)+17) = v1+2*nun;
      Val_A_rig(36*(k-1)+17) = 1;
      
      I_A_rig(36*(k-1)+18) = 6*3*(k-1)+9;
      J_A_rig(36*(k-1)+18) = v2+2*nun;
      Val_A_rig(36*(k-1)+18) = -1;
      
      b_rig(6*3*(k-1)+9) = CV(v1,3) - CV(v2,3);
      
      
      I_A_rig(36*(k-1)+19) = 6*3*(k-1)+10;
      J_A_rig(36*(k-1)+19) = w1;
      Val_A_rig(36*(k-1)+19) = 1;
      
      I_A_rig(36*(k-1)+20) = 6*3*(k-1)+10;
      J_A_rig(36*(k-1)+20) = v2;
      Val_A_rig(36*(k-1)+20) = -1;
      
      b_rig(6*3*(k-1)+10) = CV(w1,1) - CV(v2,1);
      
      I_A_rig(36*(k-1)+21) = 6*3*(k-1)+11;
      J_A_rig(36*(k-1)+21) = w1+nun;
      Val_A_rig(36*(k-1)+21) = 1;
      
      I_A_rig(36*(k-1)+22) = 6*3*(k-1)+11;
      J_A_rig(36*(k-1)+22) = v2+nun;
      Val_A_rig(36*(k-1)+22) = -1;
      
      b_rig(6*3*(k-1)+11) = CV(w1,2) - CV(v2,2);
      
      I_A_rig(36*(k-1)+23) = 6*3*(k-1)+12;
      J_A_rig(36*(k-1)+23) = w1+2*nun;
      Val_A_rig(36*(k-1)+23) = 1;
      
      I_A_rig(36*(k-1)+24) = 6*3*(k-1)+12;
      J_A_rig(36*(k-1)+24) = v2+2*nun;
      Val_A_rig(36*(k-1)+24) = -1;
      
      b_rig(6*3*(k-1)+12) = CV(w1,3) - CV(v2,3);
      
      
      I_A_rig(36*(k-1)+25) = 6*3*(k-1)+13;
      J_A_rig(36*(k-1)+25) = v1;
      Val_A_rig(36*(k-1)+25) = 1;
      
      I_A_rig(36*(k-1)+26) = 6*3*(k-1)+13;
      J_A_rig(36*(k-1)+26) = w2;
      Val_A_rig(36*(k-1)+26) = -1;
      
      b_rig(6*3*(k-1)+13) = CV(v1,1) - CV(w2,1);
      
      I_A_rig(36*(k-1)+27) = 6*3*(k-1)+14;
      J_A_rig(36*(k-1)+27) = v1+nun;
      Val_A_rig(36*(k-1)+27) = 1;
      
      I_A_rig(36*(k-1)+28) = 6*3*(k-1)+14;
      J_A_rig(36*(k-1)+28) = w2+nun;
      Val_A_rig(36*(k-1)+28) = -1;
      
      b_rig(6*3*(k-1)+14) = CV(v1,2) - CV(w2,2);
      
      I_A_rig(36*(k-1)+29) = 6*3*(k-1)+15;
      J_A_rig(36*(k-1)+29) = v1+2*nun;
      Val_A_rig(36*(k-1)+29) = 1;
      
      I_A_rig(36*(k-1)+30) = 6*3*(k-1)+15;
      J_A_rig(36*(k-1)+30) = w2+2*nun;
      Val_A_rig(36*(k-1)+30) = -1;
      
      b_rig(6*3*(k-1)+15) = CV(v1,3) - CV(w2,3);
      
      
      I_A_rig(36*(k-1)+31) = 6*3*(k-1)+16;
      J_A_rig(36*(k-1)+31) = w1;
      Val_A_rig(36*(k-1)+31) = 1;
      
      I_A_rig(36*(k-1)+32) = 6*3*(k-1)+16;
      J_A_rig(36*(k-1)+32) = w2;
      Val_A_rig(36*(k-1)+32) = -1;
      
      b_rig(6*3*(k-1)+16) = CV(w1,1) - CV(w2,1);
      
      I_A_rig(36*(k-1)+33) = 6*3*(k-1)+17;
      J_A_rig(36*(k-1)+33) = w1+nun;
      Val_A_rig(36*(k-1)+33) = 1;
      
      I_A_rig(36*(k-1)+34) = 6*3*(k-1)+17;
      J_A_rig(36*(k-1)+34) = w2+nun;
      Val_A_rig(36*(k-1)+34) = -1;
      
      b_rig(6*3*(k-1)+17) = CV(w1,2) - CV(w2,2);
      
      I_A_rig(36*(k-1)+35) = 6*3*(k-1)+18;
      J_A_rig(36*(k-1)+35) = w1+2*nun;
      Val_A_rig(36*(k-1)+35) = 1;
      
      I_A_rig(36*(k-1)+36) = 6*3*(k-1)+18;
      J_A_rig(36*(k-1)+36) = w2+2*nun;
      Val_A_rig(36*(k-1)+36) = -1;
      
      b_rig(6*3*(k-1)+18) = CV(w1,3) - CV(w2,3);

      
  end
  
  A_rig = sparse(I_A_rig,J_A_rig,Val_A_rig,6*3*edge_num_collisions,3*nun);