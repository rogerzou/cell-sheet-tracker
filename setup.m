clear all
close all
clear java

% add file paths
addpath(pwd);
addpath(genpath(fullfile(pwd, 'bin')));
addpath(genpath(fullfile(pwd, 'img')));
addpath(genpath(fullfile(pwd, 'lib')));
addpath(genpath(fullfile(pwd, 'src')));

% add dynamic java paths
javaaddpath(fullfile(pwd, 'bin'));
libpaths = genpath(fullfile(pwd, 'lib'));
allpaths = strsplit(libpaths, ':');
for ii = 1:numel(allpaths)
    javaaddpath(allpaths{ii});
end
clear
