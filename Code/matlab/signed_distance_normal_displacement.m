function [V_disp,dt,violations] = signed_distance_normal_displacement(V0,F0)
%
% Given an initial mesh (V0,F0), computes normal displacements
% V0 = V0+dt*(-N) until it finds a normal displacement that does
% not intersect the original. If it doesn't, output dt = +inf.
%
% [V_disp,dt,violations] = signed_distance_normal_displacement(V0,F0)
%
% Input:
% (V0,F0): Input mesh
% Output
% V_disp:    = V0+dt*(-N) such that (V_disp,F0) doesn't intersect (V0,F0)
% dt:     dt above
% violations:   number of times that a normal at a vertex has dot
% product with some face around it smaller than zero

% calculate normals according to signed_distance_direction function
[N_dist,~,~] = signed_distance_direction(V0,V0,F0);

% Triangle normals
N_faces = normals(V0,F0);

% compute number of violations
violations = 0;
for i=1:size(V0,1)
    faces = mod(find(F0 == i)-1,size(F0,1))+1;
    for k=1:size(faces,1)
        if (dot(N_faces(faces(k),:),N_dist(i,:))>0)
            violations = violations+1;
        end
    end
end

% Now compute normal displacements
dt = 1e-16;
IF = 1;
while (size(IF,1)>0 && dt<1)
    
    V_disp = V0 + dt*N_dist;
    IF = intersect_other(V0,F0,V_disp,F0,'FirstOnly',true);
    dt = dt*10;
    
end
% If couldn't find dt, output infinity
if (dt==1)
    V_disp = zeros(size(V0));
    dt = inf;
end

