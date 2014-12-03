function [V_coarse F_coarse] = symmetry_x_decimation(V0,F0,N)
% [V_coarse F_coarse] = symmetry_x_decimation(V0,F0,N)
% Input:
% (V0, F0): input mesh
% N: number of faces of the output mesh
% output
% (V_coarse,F_coarse): symmetric decimated mesh with N faces


% generate coarse layer with qslim and fix it
[~,V_coarse_,F_coarse_,~] = qslim(V0,F0,N);
[V_coarse_,F_coarse_] = meshfix_matlab(V_coarse_,F_coarse_);

% intersect with x=0 plane
[V_inter,F_inter,~] = selfintersect([V_coarse_; 0 -1 -1; 0 1 -1; 0 -1 1; 0 1 1],...
    [F_coarse_; size(V_coarse_,1)+1 size(V_coarse_,1)+2 size(V_coarse_,1)+3;...
    size(V_coarse_,1)+2 size(V_coarse_,1)+4 size(V_coarse_,1)+3]);

% select left side of the mesh (x<0)
b=1;
while b<=size(F_inter,1)
    if ((V_inter(F_inter(b,1),1)>=0) && (V_inter(F_inter(b,2),1)>=0) && (V_inter(F_inter(b,3),1)>=0))
        F_inter = [F_inter(1:b-1,:);F_inter(b+1:end,:)];
    else
        b = b+1;
    end
end
[RV,IM] = remove_unreferenced(V_inter,F_inter);
V_inter = RV;
F_inter = IM(F_inter);

% duplicate and glue
V_coarse_ = [V_inter; [-V_inter(:,1) V_inter(:,2) V_inter(:,3)]];
F_coarse_ = [F_inter; [size(V_inter,1)+F_inter(:,1) size(V_inter,1)+F_inter(:,3) size(V_inter,1)+F_inter(:,2)]];
[V_coarse_,~,SVJ] = remove_duplicate_vertices(V_coarse_,1e-7);
F_coarse_ = SVJ(F_coarse_);

V_coarse = V_coarse_;
F_coarse = F_coarse_;
