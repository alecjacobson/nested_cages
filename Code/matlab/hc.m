function DW = hc(DV,DT,DF,CV,CE)
  % Solve surface problem first.
  [SV,IM] = remove_unreferenced(DV,DF);
  SF = IM(DF);
  [b,bc] = boundary_conditions(SV,SF,CV,1:size(CV,1),[],CE);
  SW = kharmonic(SV,SF,b,bc);

  b = find(IM<=size(SV,1));
  bc = SW;
  DW = kharmonic(DV,DT,b,bc);
end
