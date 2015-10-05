
s = [];
for pass = 1:2
  for k = 7:-1:0
    set(gcf,'color','w');
    XV = max([ ...
      load_mesh('../../Meshes/Results/Model9_isoCGAL/Model9_7.obj'); ...
      load_mesh('../../Meshes/Results/Model9_varap/Model9_varap_7.obj')]);
    NV = min([ ...
      load_mesh('../../Meshes/Results/Model9_isoCGAL/Model9_7.obj'); ...
      load_mesh('../../Meshes/Results/Model9_varap/Model9_varap_7.obj')]);
    if pass == 1
      filename = ...
        sprintf('../../Meshes/Results/Model9_varap/Model9_varap_%d.obj',k);
    else
      filename = ...
        sprintf('../../Meshes/Results/Model9_isoCGAL/Model9_%d.obj',k);
    end
    if k==7
      t = tsurf([1 2 1],[XV;NV],fphong,'EdgeColor','none');
    end
    [V,F] = load_mesh(filename);
    [TV,TT,TF] = tetgen(V,F,'Flags','-q10');
    left_foot  = TV(:,2)<0.1 & TV(:,1)<0.49;
    right_hand = TV(:,2)<0.42 & TV(:,1)>0.675;
    R = left_foot*1 + right_hand*2;
    [b,bc] = region_boundary_conditions(R);
    L = cotmatrix3(TV,TT);
    M = massmatrix3(TV,TT,'barycentric');
    TZ = min_quad_with_fixed(-L,M*ones(size(TV,1),1),b,bc(:,1)+bc(:,2));
    Z = TZ(1:size(V,1),:);
    colormap(parula(15));
    caxis([-0.41078 1]);
    set(t,'Vertices',[V;XV;NV],'Faces',F,'FaceColor','interp','CData',mean(Z(F),2));
    if k==7
      axis equal;
      view(2);
      camproj('persp');
      set(gca,'Visible','off');
      l = light('Position',[0.2 0.8 0.5],'Style','infinite');
    end
    delete(s);
    s = add_shadow(t,l,'Ground',[0 -1 0 NV(:,2)]);
    if k==7
      axis manual;
    end
    apply_ambient_occlusion(t,'AddLights',false);
    pause();
  end
end
