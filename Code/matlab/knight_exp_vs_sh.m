%[V,F] = load_mesh('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/inward/Model4_input.obj');
%[CV,CF] = load_mesh('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/inward/Model4_coarse.obj');
%
%R = axisangle2matrix([1 0 0],-pi/2);
%NV = [];
%XV = [];
%% get bounds
%for it = 10:10:100
%  [V] =  load_mesh(sprintf('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/inward/Model4_fine_%d.obj',it));
%  [CV] = load_mesh(sprintf('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/Model4_coarse_%d.obj',it));
%  NV = min([NV;[CV;V]*R]);
%  XV = max([XV;[CV;V]*R]);
%end
%
% reset to original

for method = {'shrink','expand'}
  for it = 0:10:100
    if it == 0
      [V,F] = load_mesh('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/inward/Model4_input.obj');
      V = V*R;
      [CV,CF] = load_mesh('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/inward/Model4_coarse.obj');
      CV = CV*R;
      clf;
      hold on;
      t = tsurf(F,[V;NV;XV], ...
          'FaceVertexCData',repmat([0.3 0.4 0.8],size(F,1),1), ...
          'EdgeColor','none', ...
          'FaceAlpha',0.9,'SpecularStrength',0.0);
      tc = tsurf(CF,CV, ...
          'FaceAlpha',0.4,'FaceColor',[0.9 0.9 0.9],'LineWidth',4);
      hold off;
      view(11,8);
      camproj('persp');
      l = light('Position',[0.3 -0.3 0.8],'Style','infinite');
      s = add_shadow(t,l,'Ground',[0   0 -1 NV(:,3)]);
      sc = add_shadow(tc,l,'Ground',[0 0 -1 NV(:,3)-1e-4]);
      s.FaceColor = [1 1 1]*0.8;
      sc.FaceColor = [1 1 1]*0.95;
      axis equal;
      axis manual;
      expand_axis(1.2);
      set(gca,'Visible','off')
      set(gcf,'color','w');
      apply_ambient_occlusion(t,'AddLights',false);
    else
      switch method{1}
      case 'expand'
        CV = load_mesh(sprintf('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/Model4_coarse_%d.obj',it))*R;
        set(tc,'Vertices',CV);
      case 'shrink'
        V = load_mesh(sprintf('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/inward/Model4_fine_%d.obj',it))*R;
        set(t,'Vertices',V,'FaceVertexCData',repmat([0.3 0.4 0.8],size(F,1),1));
        apply_ambient_occlusion(t,'AddLights',false);
        delete(s);
        s = add_shadow(t,l,'Ground',[0   0 -1 NV(:,3)]);
        s.FaceColor = [1 1 1]*0.8;
      end
    end
    [~,~,IF] = selfintersect(CV,CF);
    IO = intersect_other(V,F,CV,CF);
    CMc = [nan nan nan;1.0 0.7 0.4;1.0 1.0 0.4];
    CM = [0.9 0.9 0.9;nan(2,3)];
    I = min(ismember(1:size(CF,1),IF)'+2*ismember(1:size(CF,1),IO(:,2))'+1,size(CM,1));;
    tc.FaceVertexCData = CM(I,:);
    cc = copyobj(tc,tc.Parent);
    cc.FaceVertexCData = CMc(I,:);
    cc.FaceColor = 'flat';
    cc.FaceAlpha = 0.9;
    cc.AmbientStrength = 0.4;
    tc.FaceColor = 'flat';
    delete(sc);
    sc = add_shadow(tc,l,'Ground',[0 0 -1 NV(:,3)-1e-4]);
    sc.FaceColor = [1 1 1]*0.95;
  
    drawnow;
    imwrite(myaa('raw'),sprintf('knight-%s-%03d.png',method{1},it));
    if exist('cc') && ishandle(cc) && it < 100
      delete(cc);
    end
  end
end
