function increaseVersion(bumpType)
% increaseVersion
%   Upgrades the model to a new version. Run this function after merging
%   changes to the master branch for making a new release.
%
%   bumpType    One of the following 3 strings: "major", "minor" or
%               "patch", indicating the type of increase of version to be
%               performed.
%
%   NOTE: This function requires a git wrapper added to the MATLAB search
%         path: https://github.com/manur/MATLAB-git
%
%   Usage: increaseVersion(bumpType)
%
%   Benjamin Sanchez, 2018-09-25
%

%Check if in master:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if ~strcmp(currentBranch,'master')
    error('ERROR: not in master')
end

%Bump version number:
oldModel   = load('../ModelFiles/mat/yeastGEM.mat');
oldVersion = oldModel.model.modelID;
oldVersion = oldVersion(strfind(oldVersion,'_v')+2:end);
oldVersion = str2double(strsplit(oldVersion,'.'));
newVersion = oldVersion;
switch bumpType
    case 'major'
        newVersion(1) = newVersion(1) + 1;
        newVersion(2) = 0;
        newVersion(3) = 0;
    case 'minor'
        newVersion(2) = newVersion(2) + 1;
        newVersion(3) = 0;
    case 'patch'
        newVersion(3) = newVersion(3) + 1;
    otherwise
        error('ERROR: invalid input. Use "major", "minor" or "patch"')
end
newVersion = num2str(newVersion,'%d.%d.%d');

%Check if history has been updated:
fid     = fopen('../history.md','r');
history = fscanf(fid,'%s');
fclose(fid);
if ~contains(history,['yeast' newVersion ':'])
    error('ERROR: update history.md first')
end

%Load model:
initCobraToolbox
model = readCbModel('../ModelFiles/xml/yeastGEM.xml');

%Include tag and save model:
model.modelID = ['yeastGEM_v' newVersion];
saveYeastModel(model,false)

%Check if any file changed (except for history.md and 1 line in yeastGEM.xml):
diff   = git('diff --numstat');
diff   = strsplit(diff,'\n');
change = false;
for i = 1:length(diff)
    diff_i = strsplit(diff{i},'\t');
    if length(diff_i) == 3
        %.xml file: 1 line should be added & 1 line should be deleted
        if strcmp(diff_i{3},'ModelFiles/xml/yeastGEM.xml')
            if eval([diff_i{1} ' > 1']) || eval([diff_i{2} ' > 1'])
                disp(['NOTE: File ' diff_i{3} ' is changing more than expected'])
                change = true;
            end
        %Any other file except for history.md: no changes should be detected
        elseif ~strcmp(diff_i{3},{'history.md'})
            disp(['NOTE: File ' diff_i{3} ' is changing'])
            change = true;
        end
    end
end
if change
    error(['Some files are changing from devel. To fix, first update devel, ' ...
        'then merge to master, and try again.'])
end

%Allow .mat & .xlsx storage:
copyfile('../.gitignore','backup')
fin  = fopen('backup','r');
fout = fopen('../.gitignore','w');
still_reading = true;
while still_reading
  inline = fgets(fin);
  if ~ischar(inline)
      still_reading = false;
  elseif ~startsWith(inline,'*.mat') && ~startsWith(inline,'*.xlsx')
      fwrite(fout,inline);
  end
end
fclose('all');
delete('backup');

%Store model as .mat:
save('../ModelFiles/mat/yeastGEM.mat','model');

%Convert to RAVEN format and store model as .xlsx:
model = ravenCobraWrapper(model);
model.annotation.defaultLB    = -1000;
model.annotation.defaultUB    = +1000;
model.annotation.taxonomy     = 'taxonomy/559292';
model.annotation.givenName    = 'Benjamin';
model.annotation.familyName   = 'Sanchez';
model.annotation.email        = 'bensan@chalmers.se';
model.annotation.organization = 'Chalmers University of Technology';
model.annotation.note         = 'Saccharomyces cerevisiae - strain S288C';
exportToExcelFormat(model,'../ModelFiles/xlsx/yeastGEM.xlsx');

%Update version file:
fid = fopen('../version.txt','wt');
fprintf(fid,newVersion);
fclose(fid);

end
