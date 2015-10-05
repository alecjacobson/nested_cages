[V,F] = load_mesh('/Users/ajx/Documents/nested_cages/Meshes/Results/Model4_volume/Model4_0.obj');
FV = V;
FF = F;
[FV,FF] = upsample(FV,FF);
t = tsurf(F,V,'FaceAlpha',0.2,'EdgeAlpha',0.0,'FaceColor','r');
hold on;
tf = tsurf(FF,FV,'FaceAlpha',0.2,'EdgeAlpha',0.0,'FaceColor','b');
hold off;
view(2);
axis equal;
delta = 1e-4;

S = inf;
I = 1:size(FV,1);
for iter = 1:inf
  S_prev = S;
  [D,S] = signed_distance_direction(FV(I,:),V,F);
  if abs(min(S)-min(S_prev)) < 1e-6*delta
    fprintf('Small change.\n');
    break;
  end
  P = (S<S_prev);
  if ~any(P)
    fprintf('No progress.\n');
    break;
  end
  D = D(P,:);
  S = S(P,:);
  I = I(P);
  FV(I,:) = FV(I,:) + delta*D;
  if mod(iter,20) == 1
    tf.Vertices = FV;
    drawnow;
  end
end
