function write_cages(prefix,cages_V,cages_F)
  [s,r] = system(sprintf('mkdir -p `dirname %s`',prefix));
  if s~=0
    error(r)
  end
  for m=1:numel(cages_V)-1
    writeOBJ( ...
      sprintf('%s_%d.obj',prefix,numel(cages_V)-m),cages_V{m},cages_F{m});
  end
  writeOBJ(sprintf('%s.obj',prefix),cages_V{end},cages_F{end});
end
