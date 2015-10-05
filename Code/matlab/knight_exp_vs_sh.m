[V,F] = load_mesh('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/inward/Model4_input.obj');
[CV,CF] = load_mesh('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/inward/Model4_coarse.obj');

R = axisangle2matrix([1 0 0],-pi/2);
NV = [];
XV = [];
% get bounds
for it = 10:10:100
  [V] =  load_mesh(sprintf('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/inward/Model4_fine_%d.obj',it));
  [CV] = load_mesh(sprintf('/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4/Model4_coarse_%d.obj',it));
  NV = min([NV;[CV;V]*R]);
  XV = max([XV;[CV;V]*R]);
end

% reset to original

prefix = '/Users/ajx/Documents/nested_cages/Meshes/Expansion/Model4';
data = [];
for method = {'shrink','expand'}
  for it = 0:10:100
    if it == 0
      [V,F] = load_mesh([prefix '/inward/Model4_input.obj']);
      V = V*R;
      [CV,CF] = load_mesh([prefix '/inward/Model4_coarse.obj']);
      CV = CV*R;
      clf;
      data = [];
    else
      switch method{1}
      case 'expand'
        CV = load_mesh(sprintf('%s/Model4_coarse_%d.obj',prefix,it))*R;
      case 'shrink'
        V = load_mesh(sprintf('%s/inward/Model4_fine_%d.obj',prefix,it))*R;
      end
    end
    data = render_in_cage(V,F,CV,CF, ...
      'BoundingBox',[NV;XV],'Data',data,'AmbientOcclusion',true);

    imwrite(myaa('raw'),sprintf('knight-%s-%03d.png',method{1},it));
  end
end
