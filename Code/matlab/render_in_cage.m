function data = render_in_cage(V,F,CV,CF,varargin)

  % Delete a handle, but only if it's not empty and actually a handle
  %
  % Inputs:
  %   h  handle id
  function guarded_delete(h)
    if ~isempty(h) && ishandle(h)
      delete(h)
    end
  end

  blue = [0.0482,0.3651,0.8722];
  vaa = [11 8];
  l_position = [0.3 -0.3 0.8];
  ground = [];
  BB = [];
  ao = true;
  color_intersections = true;
  data = [];
  ao_factor = 1.0;

  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    { 'AmbientOcclusion', ...
      'AOFactor', ...
      'BoundingBox', ...
      'ColorIntersections', ...
      'Ground', ...
      'LightPosition', ...
      'View', ...
      'Data' ...
    }, ...
    { 'ao', ...
      'ao_factor', ...
      'BB', ...
      'color_intersections', ...
      'ground', ...
      'l_position', ...
      'vaa', ...
      'data' ...
    });
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  if isempty(data)
    data.t = [];
    data.l = [];
    data.tc = [];
    data.s = [];
    data.sc = [];
    data.cc = [];
    data.A = [];
  end

  if isempty(BB)
    BB = [min([CV;V]);max([CV;V])];
  end

  if isempty(ground)
    ground = [0   0 -1 BB(1,3)];
  end

  if isempty(data.t)
    data.t = tsurf(F,[V;BB], ...
        'FaceVertexCData',repmat(blue,size(F,1),1), ...
        'EdgeColor','none', ...
        'FaceAlpha',0.9,'SpecularStrength',0.0);
  else
    if isempty(V)
      V = data.t.Vertices;
    else
      set(data.t,'Vertices',V,'FaceVertexCData',repmat(blue,size(F,1),1));
    end
    if isempty(F)
      F = data.t.Faces;
    else
      set(data.t,'Faces',F);
    end
  end
  old_hold = ishold;
  hold on;
  if isempty(data.tc)
    data.tc = tsurf(CF,CV, ...
        'FaceAlpha',0.4,'FaceColor',[0.9 0.9 0.9],'LineWidth',4);
  else
    if isempty(CV)
      CV = data.tc.Vertices;
    else
      set(data.tc,'Vertices',CV);
    end
    if isempty(CF)
      CF = data.tc.Faces;
    else
      set(data.tc,'Faces',CF);
    end
  end
  hold([repmat('on',old_hold,old_hold) repmat('off',~old_hold,~old_hold)])
  
  view(vaa);
  camproj('persp');
  guarded_delete(data.l);
  data.l = light('Position',l_position,'Style','infinite');
  guarded_delete(data.s);
  data.s  = add_shadow(data.t,data.l, 'Ground',ground);
  data.s.FaceColor = [1 1 1]*0.8;


  if color_intersections
    [~,~,IF] = selfintersect(CV,CF);
    IO = intersect_other(V,F,CV,CF);
    CMc = [nan nan nan;1.0 0.7 0.4;1.0 1.0 0.4];
    CM = [0.9 0.9 0.9;nan(2,3)];
    I = min(ismember(1:size(CF,1),IF)'+2*ismember(1:size(CF,1),IO(:,2))'+1,size(CM,1));;
    data.tc.FaceVertexCData = CM(I,:);
    guarded_delete(data.cc);
    data.cc = copyobj(data.tc,data.tc.Parent);
    data.cc.FaceVertexCData = CMc(I,:);
    data.cc.FaceColor = 'flat';
    data.cc.FaceAlpha = 0.9;
    data.cc.AmbientStrength = 0.4;
    data.tc.FaceColor = 'flat';
  end

  guarded_delete(data.sc);
  data.sc = add_shadow(data.tc,data.l,'Ground',ground-[0 0 0 1e-3]);
  data.sc.FaceColor = [1 1 1]*0.95;

  if isempty(data.A)
    axis equal;
    axis manual;
    expand_axis(1.2);
    data.A = axis;
  else
    axis(data.A);
  end
  set(gca,'Visible','off')
  set(gcf,'color','w');

  if ao
    apply_ambient_occlusion(data.t,'AddLights',false,'Factor',ao_factor);
  end

end
