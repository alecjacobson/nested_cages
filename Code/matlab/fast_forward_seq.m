prefix = '../../FastForward/MeshSequences/fertility/';
% retrive number of meshes
n = 0;

while true
  if exist(sprintf('%s/V_%d.obj',prefix,n+1),'file') ~= 2
    break;
  end
  n = n+1;
end

% global rotation
R = axisangle2matrix([1 0 0],-pi/2);
% retrive bounds
[CV,CF] = load_mesh(sprintf('%s/CV_%d.obj',prefix,n));
CV = CV*R;
BB = [min(CV)-1e-4;max(CV)+1e-4];

suffix = @(n) [repmat('input',n==0,n==0) repmat(num2str(n),n~=0,n~=0)];

clf;drawnow;
T = get(gca,'tightinset');set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);
data = [];
for direction = {'shrink','expand'}
  for iter = 0:n
    switch direction{1}
    case 'shrink'
      [CV,CF] = load_mesh(sprintf('%s/CV_%s.obj',prefix,suffix(0)));
      [V,F] = load_mesh(sprintf('%s/V_%s.obj',prefix,suffix(iter)));
      titer = iter;
    case 'expand'
      [CV,CF] = load_mesh(sprintf('%s/CV_%s.obj',prefix,suffix(n-iter)));
      [V,F] = load_mesh(sprintf('%s/V_%s.obj',prefix,suffix(n-iter)));
      titer = n+iter;
    end
    V = V*R;
    CV = CV*R;
    data = render_in_cage(V,F,CV,CF, ...
      'BoundingBox',BB, ...
      'Data',data,'View',[26 24],'LightPosition',[0.1 -0.1 0.8], ...
      'AOFactor',1.0, ...
      'AmbientOcclusion',true);
    drawnow;
    frame = getframe(gcf);
    imwrite(frame.cdata,sprintf('fertility-%03d.png',titer));
  end
end
