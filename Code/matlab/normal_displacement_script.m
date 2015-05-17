function no_displace = normal_displacement_script(basename)

no_displace = '';

allFiles = dir(basename);
names = {allFiles.name};
num_files = size(names,2);
for k=1:num_files
    filename = names(k);
    off_cmp = strfind(filename,'.off');
    obj_cmp = strfind(filename,'.obj');
    if (~isempty(off_cmp{1})||~isempty(obj_cmp{1}))
        fullname = sprintf('%s/%s',basename,filename{1})
        [V0,F0] = load_mesh(fullname);
        [V_disp,dt,violations] = signed_distance_normal_displacement(V0,F0);
        if dt==Inf
%             no_displace = [no_displace; fullname];
            no_displace = sprintf('%s \n %s',no_displace,fullname)
        end
    end
end
            