%levels = [1 16 26 34];
%prefix = '/Users/ajx/Documents/nested_cages/Meshes/Results/flow_comparison_2/ours_disppath_all/';
%scheme = 'ours'
levels = [16 26 36 46 56 100];
prefix = '/Users/ajx/Documents/nested_cages/Meshes/Results/flow_comparison_2/cmcf_disppath_all/';
scheme = 'cmcf';
clf;
hold on;
 t = tsurf([1 1 1],[0 0 0],'FaceColor',[0.0482,0.3651,0.8722],'EdgeColor','none');
tc = tsurf([1 1 1],[0 0 0],'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.2,'EdgeAlpha',0.5,'LineWidth',2);
l = light('Position',[0.3 -0.4 0.8],'Style','infinite');
hold off;
set(gca,'Visible','off');
set(gcf,'Color','w');
axis equal;
camproj('persp');
R = axisangle2matrix([1 0 0],-pi/2);
s = [];
for iter = levels 
  for pass = 1:2
    switch pass
    case 1
      [CV,CF] = load_mesh(sprintf('%sCV_%d.obj',prefix,iter));
    case 2
      [CV,CF] = load_mesh(sprintf('%sCV_0.obj',prefix));
    end
    [V,F] = load_mesh(sprintf('%sV_%d.obj',prefix,iter));
    [V,F] = load_mesh('/Users/ajx/Documents/nested_cages/Meshes/Results/rampant_varap_new/rampant_0.obj');
    set( t,'Vertices', V*R,'Faces',F);
    set(tc,'Vertices',CV*R,'Faces',CF);
    delete(s);
    ground = [0 0 -1 -0.3];
    s = add_shadow(t,l,'Ground',ground);
    sc = add_shadow(tc,l,'Ground',ground-[0 0 0 1e-4]);
    s.FaceColor = [1 1 1]*0.8;
    sc.FaceColor = [1 1 1]*0.95;
    if iter == levels(1)  && pass == 2
      axis manual;
      view(-21,8);
    end
    if pass == 2
    apply_ambient_occlusion(t,'AddLights',false);
    %imwrite(myaa('raw'),sprintf('rampant-alt-%s-%d-%d.png',scheme,iter,pass));
    imwrite(myaa('raw'),'rampant-alt-input.png');
    error
    end
    %pause
  end
end
