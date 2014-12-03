function EE = CTCD_edge_edge(V0,V1,F)
  % CTCD_EDGE_EDGE
  % CTCD_edge_edge(V0,V1,F)
  %
  % Continuous time collision detection (edge-edge) only
  %
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions at the begining of
  %   the time step
  %   V1  (#vertices)x3 list of mesh vertex positions at the end of
  %   the time step
  %   F  (#faces)x3 vertex indices for each face of both meshes
  %     
  % Output:
  % EE: one row per collision between two edges
  %     format: e1_v1 e1_v2 e2_v1 e2_v2 t
  %     e1_v1, e1_v2: vertices ids of edge 1, 
  %     e2_v1, e2_v2: vertices ids of edge 2, 
  %     t: time [0-1]
  
  E = edges(F);
  EE = [];
  
  tic
  for i = 1:size(E,1)
      i
      for j = i:size(E,1)
          
          q0 = V0(E(i,1),:);
          p0 = V0(E(i,2),:);
          q1 = V0(E(j,1),:);
          p1 = V0(E(j,2),:);
          
          q0end = V1(E(i,1),:);
          p0end = V1(E(i,2),:);
          q1end = V1(E(j,1),:);
          p1end = V1(E(j,2),:);
          
%           % Tyson's code
%           exact_ccd_mex(q0,p0,q1,p1,q0end,p0end,q1end,p1end,1);
          
          % Etienne's code
          t = CTCD_mex(q0,p0,q1,p1,q0end,p0end,q1end,p1end);
          
          if (t>0 && t<=1)
              EE = [EE; E(i,1) E(i,2) E(j,1) E(j,2) t]
          end
          
      end
      
  end
  toc