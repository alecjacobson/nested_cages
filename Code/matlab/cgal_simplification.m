function [W,G] = cgal_simplification(V,F,t)
  % Very simple wrapper for CGAL simplification (edge-colapse)
  %
  % [W,G] = cgal_simplification(V,F,t,'ParameterName',ParameterValue)
  %
  % Inputs:
  %   V  #V by 3 input mesh vertex positions
  %   F  #F by 3 input mesh triangle indices (1-indexed)
  %   t  target number of faces {0 for all collapse}
  % Outputs:
  %   W  #W by 3 collapsed mesh vertex positions
  %   G  #G by 3 collapsed mesh triangle indices (1-indexed)

  prefix = tempprefix();
  input = [prefix '_in.off'];
  output = [prefix '_out.off'];

  % Write input to file
  writeOFF(input,V,F);

    % prepare command string
    command = sprintf('%s %s %d %s', ...
      '/usr/bin/edge_collapse_enriched_polyhedron ',...
      input,t/size(F,1),output);
    [status,result] = system(command);
    if status ~= 0
      error(result);
    end

    [W,G] = readOFF(output);

    delete(input);
    delete(output);


end
