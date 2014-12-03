[V,F] = load_mesh('/Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/examples/shared/decimated-knight.obj');
[U,G] = load_mesh('/Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/examples/shared/cheburashka.off');
W = U;

tsurf(F,V,'FaceAlpha',0.25,'FaceColor','r');
hold on;
t = tsurf(G,U,'FaceAlpha',0.125,'FaceColor','b');
D = 0*U;
q = quiver3(U(:,1),U(:,2),U(:,3),D(:,1),D(:,2),D(:,3),0);
hold off;
axis equal;
view(2);

[sqrD_prev] = point_mesh_squared_distance(W,V,F);
w = winding_number(V,F,W);
sqrD_prev = (1-2*w).*sqrD_prev; 

step = 0.05*ones(size(U,1),1);
sqrD = sqrD_prev;
while true
  W_prev = W;
  sqrD_prev = sqrD;
  % binary search on step
  D = signed_distance_direction(W,V,F);
  while true
    D_step = bsxfun(@times,step,D);
    W = W_prev+D_step;
    set(q, ... 
      'XData',W_prev(:,1), 'YData',W_prev(:,2), 'ZData',W_prev(:,3), ...
      'UData',D_step(:,1), 'VData',D_step(:,2), 'WData',D_step(:,3));
    sqrD = point_mesh_squared_distance(W,V,F);
    w = winding_number(V,F,W);
    sqrD = (1-2*w).*sqrD; 

    too_far = sqrD>sqrD_prev;
    if any(too_far)
      step(too_far) = 0.5*step(too_far);
    else
      break;
    end
    drawnow;
    break;
  end
  drawnow;
  input('');
  set(t,'Vertices',W);

end
