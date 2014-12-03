function plot_struct = flow_plot_control(plot_struct,vis)
  % FLOW_PLOT_CONTROL
  % 
  % plot_struct = flow_plot_control(plot_struct,vis)
  %
  % Used to control the flow plot from shrink_fine_to_inside_coarse_3D.
  % Allows options such as displaying/removing quadrature points, gradients at
  % quadrature points, gradients at vertices, intersections between fine
  % and coarse mesh, etc.
  %
  % Input: 
  %   plot_struct: structure that encapsulates all data needed for plotting
  %         .V_shrink: shrunk mesh vertices
  %         .F_shrink: shrunk mesh faces
  %         .V_to_intersect: all coarse layers vertices
  %         .F_to_intersect: all coarse layers faces
  %         .V_coarse: finest coarse layers (the ones that indeces the
  %         sugned distance field) vertices
  %         .F_coarse: finest coarse layers (the ones that indeces the
  %         sugned distance field) faces
  %         .quad_points: quadrature points on the fine mesh triangles
  %         .grad_vertices: gradient at fine mesh vertices
  %         .grad_quad: gradient at quadrature points
  %         .quad_closest: quadrature points closest points on the coarse
  %         mesh
  %         .V_shrink_prev: shrunk mesh vertices at the previous step
  %         .IF: intersecting triangles (first columns are indices of
  %         coarse mesh triangles and second column of fine mesh triangles)
  %         .pc: plot handle for the coarse mesh
  %         .pc1: plot handle for the finest coarse mesh
  %         .pc_alpha: alpha channel for the coarse mesh
  %         .pv: plot handle for the fine mesh
  %         .pv_alpha: alpha channel for the fine mesh
  %         .p_quiver: quiver handle for gradients at vertices of the fine
  %         mesh
  %         .p_quiver_quad: quiver handle for gradients at quadrature
  %         points
  %         .p_close: plot handles for quadrature points closest points on
  %         the coarse mesh
  %         .p_quad: plot handle for quadrature points.
  %  vis: either true (visualization is on) or false (vis is off)
  % Output:
  %   plot_struct: modified plot_struct
  
%     plot_struct = struct('V_shrink',V_shrink,'F_shrink',F_shrink,'V_to_intersect',V_to_intersect,...
%       'F_to_intersect',F_to_intersect,'V_coarse',V_coarse,'F_coarse',F_coarse, 'quad_points',[],...
%       'grad_vertices',[],'grad_quad',[],'quad_closest',[],'V_shrink_prev',V0,'IF',[],'pc',pc,'pc1',[],'pc_alpha',0.05,'pv',pv,'pv_alpha',0.4,...
%       'p_quiver',[],'p_quiver_quad',[],'p_close',[],'p_quad',[]);
  
  hold on
  
  % plot gradients and closest points on the mesh from the previous step
  if ~isempty(plot_struct.p_quiver)
      delete(plot_struct.p_quiver)
      plot_struct.p_quiver = quiver3(plot_struct.V_shrink_prev(:,1),plot_struct.V_shrink_prev(:,2),plot_struct.V_shrink_prev(:,3),...
          10*plot_struct.grad_vertices(:,1),10*plot_struct.grad_vertices(:,2),10*plot_struct.grad_vertices(:,3),0);
      set(plot_struct.p_quiver,'Color',[0 0.75 1.0],'Linewidth',1.5);
  end
  if ~isempty(plot_struct.p_quiver_quad)
      delete(plot_struct.p_quiver_quad)
      plot_struct.p_quiver_quad = quiver3(plot_struct.quad_points(:,1),plot_struct.quad_points(:,2),plot_struct.quad_points(:,3),...
          0.01*plot_struct.grad_quad(:,1),0.01*plot_struct.grad_quad(:,2),0.01*plot_struct.grad_quad(:,3),0);
      set(plot_struct.p_quiver_quad,'Color','m','Linewidth',1.5);
  end
  if ~isempty(plot_struct.p_close)
      delete(plot_struct.p_close)
      plot_struct.p_close = plot3(plot_struct.quad_closest(:,1),plot_struct.quad_closest(:,2),...
                      plot_struct.quad_closest(:,3),'g.','markersize',10);
  end  
  if ~isempty(plot_struct.p_quad)
      delete(plot_struct.p_quad);
      plot_struct.p_quad = plot3(plot_struct.quad_points(:,1),plot_struct.quad_points(:,2),plot_struct.quad_points(:,3),'r.','markersize',10);
  end
  
  % display options and get input from keyboard
  if (~vis)
	  str_input = '';
  else
      str_input = input('[c]: all coarse meshes \n[c1]: only finest coarse level \n[f]: fine mesh \n[q]: quadrature points \n[g]: gradient at vertices \n[gq]: gradient at quadrature points \n[t]: closest points on the coarse mesh to quadrature points \n[i]: show intersecting triangles in different colors \n','s');
  end
    
  % while some string is entered, process it. 
  % If 'return' is typed with no string, moves to the next step of the flow
  while(~strcmp(str_input,''))
      
      % string starts with 'c', so we're gonna change the coarse mesh plot
      if (strcmp(str_input(1),'c'))
          % if only 'c' was entered, it is just to hide/display
          % the coarse mesh
          if (size(str_input,2)==1)
              if isempty(plot_struct.pc)
                  delete(plot_struct.pc);
                  plot_struct.pc = trisurf(plot_struct.F_to_intersect,plot_struct.V_to_intersect(:,1),...
                      plot_struct.V_to_intersect(:,2),plot_struct.V_to_intersect(:,3),...
                      'FaceColor',[0.5 0.0 0.0],'FaceAlpha',plot_struct.pc_alpha);
              else
                  delete(plot_struct.pc);
                  plot_struct.pc = [];
              end
          % elseif followed by a '=', get value of the rest of the string
          % for the alpha channel of the coarse mesh
          elseif (strcmp(str_input(2),'='))
              % if it's in good alhpa range, plot mesh with new aplha
              if (str2double(str_input(3:end))>=0.0 && str2double(str_input(3:end))<=1.0)
                  plot_struct.pc_alpha = str2double(str_input(3:end));
                  if (~isempty(plot_struct.pc))
                      delete(plot_struct.pc);
                      plot_struct.pc = trisurf(plot_struct.F_to_intersect,plot_struct.V_to_intersect(:,1),...
                          plot_struct.V_to_intersect(:,2),plot_struct.V_to_intersect(:,3),...
                          'FaceColor',[0.5 0.0 0.0],'FaceAlpha',plot_struct.pc_alpha);
                  elseif (~isempty(plot_struct.pc1))
                      delete(plot_struct.pc1);
                      plot_struct.pc1 = trisurf(plot_struct.F_coarse,plot_struct.V_coarse(:,1),...
                        plot_struct.V_coarse(:,2),plot_struct.V_coarse(:,3),...
                        'FaceColor',[0.5 0.0 0.0],'FaceAlpha',plot_struct.pc_alpha);
                  end
              % else display a message and do nothing
              else
                  disp('please specify a proper alpha value for the coarse mesh');
              end
          % elseif followed by a '1', display only the finest coarse mesh
          % (the one that induces the signed distance field)
          elseif (strcmp(str_input(2),'1'))
              if isempty(plot_struct.pc1)
                  delete(plot_struct.pc1);
                  plot_struct.pc1 = trisurf(plot_struct.F_coarse,plot_struct.V_coarse(:,1),...
                      plot_struct.V_coarse(:,2),plot_struct.V_coarse(:,3),...
                      'FaceColor',[0.5 0.0 0.0],'FaceAlpha',plot_struct.pc_alpha);
                  delete(plot_struct.pc);
                  plot_struct.pc = [];
              else
                  delete(plot_struct.pc);
                  plot_struct.pc = trisurf(plot_struct.F_to_intersect,plot_struct.V_to_intersect(:,1),...
                      plot_struct.V_to_intersect(:,2),plot_struct.V_to_intersect(:,3),...
                      'FaceColor',[0.5 0.0 0.0],'FaceAlpha',plot_struct.pc_alpha);
                  delete(plot_struct.pc1);
                  plot_struct.pc1 = [];
              end
          end
          
      % string starts with 'f', so we're gonna change the fine mesh plot
      elseif (strcmp(str_input(1),'f'))
          % if only 'f' was entered, it is just to hide/display
          % the fine mesh
          if (size(str_input,2)==1)
              if isempty(plot_struct.pv)
                  delete(plot_struct.pv);
                  plot_struct.pv = trisurf(plot_struct.F_shrink,plot_struct.V_shrink(:,1),...
                    plot_struct.V_shrink(:,2),plot_struct.V_shrink(:,3),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',plot_struct.pv_alpha);
              else
                  delete(plot_struct.pv);
                  plot_struct.pv = [];
              end
          % elseif followed by a '=', get value of the rest of the string
          % for the alpha channel of the fine mesh
          elseif (strcmp(str_input(2),'='))
              % if it's in good alhpa range, plot mesh with new aplha
              if (str2double(str_input(3:end))>=0.0 && str2double(str_input(3:end))<=1.0)
                  plot_struct.pv_alpha = str2double(str_input(3:end));
                  delete(plot_struct.pv);
                  plot_struct.pv = trisurf(plot_struct.F_shrink,plot_struct.V_shrink(:,1),...
                        plot_struct.V_shrink(:,2),plot_struct.V_shrink(:,3),...
                        'FaceColor',[0.0 0.0 0.8],'FaceAlpha',plot_struct.pv_alpha);
              % else display a message and do nothing
              else
                  disp('please specify a proper alpha value for the coarse mesh');
              end
          elseif (strcmp(str_input(2),'w'))
              if (plot_struct.show_weights==0)
                  plot_struct.show_weights = 1;
              elseif (plot_struct.show_weights==1)
                  plot_struct.show_weights = 0;
              end
          end
          
      % string starts with 'q', so we're gonna hide/display quadrature points
      elseif (strcmp(str_input(1),'q'))
          % if only 'q' was entered, it is just to remove/display
          % the quadrature points
          if (size(str_input,2)==1)
              if isempty(plot_struct.p_quad)
                  delete(plot_struct.p_quad);
                  plot_struct.p_quad = plot3(plot_struct.quad_points(:,1),plot_struct.quad_points(:,2),...
                      plot_struct.quad_points(:,3),'r.','markersize',10);
              else
                  delete(plot_struct.p_quad);
                  plot_struct.p_quad = [];
              end
          end
          
      % 'g' means we're gonna hids/display gradients at vertices, 
      % 'gq' gradients at quadrature points
      elseif (strcmp(str_input(1),'g'))
          % if only 'g' was entered, it is just to hide/display
          % the gradients
          if (size(str_input,2)==1)
              if isempty(plot_struct.p_quiver)
                  disp('\n \n ATTENTION: displaying 10*(steps.*grad_V) \n \n')
                  delete(plot_struct.p_quiver);
                  plot_struct.p_quiver = quiver3(plot_struct.V_shrink_prev(:,1),plot_struct.V_shrink_prev(:,2),plot_struct.V_shrink_prev(:,3),...
                        10*plot_struct.grad_vertices(:,1),10*plot_struct.grad_vertices(:,2),10*plot_struct.grad_vertices(:,3),0);
                  set(plot_struct.p_quiver,'Color',[0 0.75 1.0],'Linewidth',1.5);
              else
                  delete(plot_struct.p_quiver);
                  plot_struct.p_quiver = [];
              end
          elseif (strcmp(str_input(2),'q'))
              if isempty(plot_struct.p_quiver_quad)
                  disp('\n \n ATTENTION: displaying 0.01*(steps.*grad_V) \n \n')
                  delete(plot_struct.p_quiver_quad);
                  plot_struct.p_quiver_quad = quiver3(plot_struct.quad_points(:,1),plot_struct.quad_points(:,2),plot_struct.quad_points(:,3),...
                        0.01*plot_struct.grad_quad(:,1),0.01*plot_struct.grad_quad(:,2),0.01*plot_struct.grad_quad(:,3),0);
                  set(plot_struct.p_quiver_quad,'Color','m','Linewidth',1.5);
              else
                  delete(plot_struct.p_quiver_quad);
                  plot_struct.p_quiver_quad = [];
              end
          end
          
      % 't' means we're gonna display/hide closest points on the coarse mesh    
      elseif (strcmp(str_input(1),'t'))
          if (size(str_input,2)==1)
              if isempty(plot_struct.p_close)
                  delete(plot_struct.p_close);
                  plot_struct.p_close = plot3(plot_struct.quad_closest(:,1),plot_struct.quad_closest(:,2),...
                      plot_struct.quad_closest(:,3),'g.','markersize',10);
              else
                  delete(plot_struct.p_close);
                  plot_struct.p_close = [];
              end
          end
          
      % 'i' means we're gonna display/hide intersections between meshes   
      elseif (strcmp(str_input(1),'i'))
          if (size(str_input,2)==1)
              % if there is no intersection data, calculate the
              % intersections and plot them
              if (size(plot_struct.IF,1)==0)
                  plot_struct.IF = intersect_other(plot_struct.V_to_intersect,plot_struct.F_to_intersect,...
                      plot_struct.V_shrink_prev,plot_struct.F_shrink);
                  int_coarse = zeros(size(plot_struct.F_to_intersect,1),1);
                  int_coarse(plot_struct.IF(:,1)) = 1;
                  int_fine = zeros(size(plot_struct.F_shrink,1),1);
                  int_fine(plot_struct.IF(:,2)) = 1;
                  if (~isempty(plot_struct.pc))
                      delete(plot_struct.pc);
                      plot_struct.pc = trisurf(plot_struct.F_to_intersect,plot_struct.V_to_intersect(:,1),...
                          plot_struct.V_to_intersect(:,2),plot_struct.V_to_intersect(:,3),...
                          int_coarse,'FaceAlpha',plot_struct.pc_alpha);
                  end
                  if (~isempty(plot_struct.pc1))
                      delete(plot_struct.pc1);
                      plot_struct.pc1 = trisurf(plot_struct.F_coarse,plot_struct.V_coarse(:,1),...
                      plot_struct.V_coarse(:,2),plot_struct.V_coarse(:,3),...
                      int_coarse(1:size(plot_struct.F_coarse,1)),'FaceAlpha',plot_struct.pc_alpha);
                  end
                  if (~isempty(plot_struct.pv))
                      delete(plot_struct.pv);
                      plot_struct.pv = trisurf(plot_struct.F_shrink,plot_struct.V_shrink_prev(:,1),...
                          plot_struct.V_shrink_prev(:,2),plot_struct.V_shrink_prev(:,3),...
                          int_fine,'FaceAlpha',plot_struct.pv_alpha);
                  end
              % if there is intersection data, remove it
              else
                  
                  plot_struct.IF = [];
                  
                  if (~isempty(plot_struct.pc))
                      delete(plot_struct.pc);
                      plot_struct.pc = trisurf(plot_struct.F_to_intersect,plot_struct.V_to_intersect(:,1),...
                          plot_struct.V_to_intersect(:,2),plot_struct.V_to_intersect(:,3),...
                          'FaceColor',[0.5 0.0 0.0],'FaceAlpha',plot_struct.pc_alpha);
                  end
                  if (~isempty(plot_struct.pc1))
                      delete(plot_struct.pc1);
                      plot_struct.pc1 = trisurf(plot_struct.F_coarse,plot_struct.V_coarse(:,1),...
                        plot_struct.V_coarse(:,2),plot_struct.V_coarse(:,3),...
                        'FaceColor',[0.0 0.0 0.8],'FaceAlpha',plot_struct.pc_alpha);
                  end
                  if (~isempty(plot_struct.pv))
                      delete(plot_struct.pv);
                      plot_struct.pv = trisurf(plot_struct.F_shrink,plot_struct.V_shrink_prev(:,1),...
                          plot_struct.V_shrink_prev(:,2),plot_struct.V_shrink_prev(:,3),...
                          'FaceColor',[0.0 0.0 0.8],'FaceAlpha',plot_struct.pv_alpha);
                  end
                  
              end

          end
              
              
      end
      
      % wait for another string or 'return' to move to the next flow step
      str_input = input('[c]: coarse mesh \n[c1]: only finest coarse level \n[f]: fine mesh \n[q]: quadrature points \n[g]: gradient at vertices \n[gq]: gradient at quadrature points \n[t]: closest points on the coarse mesh to quadrature points \n[i]: show intersecting triangles in different colors \n','s');
      
  end
  
  % plot new mesh with plot settings
  
  % plot fine mesh with intersection data depending if it set
  if (~isempty(plot_struct.pv)&&(size(plot_struct.IF,1)==0))
      delete(plot_struct.pv);
      plot_struct.pv = trisurf(plot_struct.F_shrink,plot_struct.V_shrink(:,1),...
          plot_struct.V_shrink(:,2),plot_struct.V_shrink(:,3),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',plot_struct.pv_alpha);
  elseif (~isempty(plot_struct.pv)&&(size(plot_struct.IF,1)>0))
      delete(plot_struct.pv);
      plot_struct.IF = intersect_other(plot_struct.V_to_intersect,plot_struct.F_to_intersect,...
          plot_struct.V_shrink,plot_struct.F_shrink);
      int_fine = zeros(size(plot_struct.F_shrink,1),1);
      int_fine(plot_struct.IF(:,2)) = 1;
      plot_struct.pv = trisurf(plot_struct.F_shrink,plot_struct.V_shrink(:,1),...
          plot_struct.V_shrink(:,2),plot_struct.V_shrink(:,3),...
          int_fine,'FaceAlpha',plot_struct.pv_alpha);
  end
  % plot coarse mesh with intersection data depending if it set
  if (~isempty(plot_struct.pc)&&(size(plot_struct.IF,1)==0))
      delete(plot_struct.pc);
      plot_struct.pc = trisurf(plot_struct.F_to_intersect,plot_struct.V_to_intersect(:,1),...
          plot_struct.V_to_intersect(:,2),plot_struct.V_to_intersect(:,3),'FaceColor',[0.5 0.0 0.0],'FaceAlpha',plot_struct.pc_alpha);
  elseif (~isempty(plot_struct.pc)&&(size(plot_struct.IF,1)>0))
      delete(plot_struct.pc);
      plot_struct.IF = intersect_other(plot_struct.V_to_intersect,plot_struct.F_to_intersect,...
          plot_struct.V_shrink,plot_struct.F_shrink);
      int_coarse = zeros(size(plot_struct.F_to_intersect,1),1);
      int_coarse(plot_struct.IF(:,1)) = 1;
      plot_struct.pc = trisurf(plot_struct.F_to_intersect,plot_struct.V_to_intersect(:,1),...
             plot_struct.V_to_intersect(:,2),plot_struct.V_to_intersect(:,3),...
             int_coarse,'FaceAlpha',plot_struct.pc_alpha);
      
  end
  % plot finest coarse mesh with intersection data depending if it set
  if (~isempty(plot_struct.pc1)&&(size(plot_struct.IF,1)==0))
      delete(plot_struct.pc1);
      plot_struct.pc1 = trisurf(plot_struct.F_coarse,plot_struct.V_coarse(:,1),...
          plot_struct.V_coarse(:,2),plot_struct.V_coarse(:,3),'FaceColor',[0.5 0.0 0.0],'FaceAlpha',plot_struct.pc_alpha);
  elseif (~isempty(plot_struct.pc1)&&(size(plot_struct.IF,1)>0))
      delete(plot_struct.pc1);
      plot_struct.IF = intersect_other(plot_struct.V_to_intersect,plot_struct.F_to_intersect,...
          plot_struct.V_shrink,plot_struct.F_shrink);
      int_coarse = zeros(size(plot_struct.F_to_intersect,1),1);
      int_coarse(plot_struct.IF(:,1)) = 1;
      plot_struct.pc1 = trisurf(plot_struct.F_coarse,plot_struct.V_coarse(:,1),...
           plot_struct.V_coarse(:,2),plot_struct.V_coarse(:,3),...
           int_coarse(1:size(plot_struct.F_coarse,1)),'FaceAlpha',plot_struct.pc_alpha);
  end
  
  drawnow;
  hold off;