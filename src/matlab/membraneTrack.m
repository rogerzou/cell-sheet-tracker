function [data] = membraneTrack(ALL, GT, options)
%MEMBRANETRACK Track graph structure in image stream.
%
% INPUTS
% ALL: 3D image stack.
% GT: ground truth.
% options: options struct.
%
% CALLEE functions
%   embryoInitGraph
%   updateGraph
%   optIterGraph
%
% @author Roger Zou
% @date 8/15/15

N = size(ALL,3);
l = 17; w = 25; alpha = 1; interval = 0; spacing = 20;
verboseE = false; verboseG = false; siftflow = false; parallel = false;
edgetype = 'A';
fname = ['tmp_', datestr(clock)];
if nargin==nargin('membraneTrack') && ~isempty(options)
    if any(strcmp('l',fieldnames(options)))
        l = options.l;
    end
    if any(strcmp('w',fieldnames(options)))
        w = options.w;
    end
    if any(strcmp('alpha',fieldnames(options)))
        alpha = options.alpha;
    end
    if any(strcmp('interval',fieldnames(options)))
        interval = options.interval;
    end
    if any(strcmp('spacing',fieldnames(options)))
        spacing = options.spacing;
    end
    if any(strcmp('verboseE',fieldnames(options)))
        verboseE = options.verboseE;
    end
    if any(strcmp('verboseG',fieldnames(options)))
        verboseG = options.verboseG;
    end
    if any(strcmp('parallel',fieldnames(options)))
        parallel = options.parallel;
    end
    if any(strcmp('siftflow',fieldnames(options)))
        siftflow = options.siftflow;
    end
    if any(strcmp('edgetype',fieldnames(options)))
        edgetype = options.edgetype;
    end
    if any(strcmp('fname',fieldnames(options)))
        fname = options.fname;
    end
end
if edgetype=='C' || edgetype=='D'
    ALL = max(ALL(:)) - ALL;
end

% optionally start parallel pool
if parallel && isempty(gcp('nocreate'))
    parpool
end

% set up some optimization structs
optOptions = struct('parallel',parallel,'verboseE',verboseE,'verboseG',verboseG,'siftflow',siftflow);
structC = struct('alpha',alpha,'l',l,'w', w,'interval',interval,'spacing',spacing);

% select edge type
if edgetype=='A'
    structD = getStructD_A;
elseif edgetype=='B'
    structD = getStructD_B;
elseif edgetype=='C'
    structD = getStructD_C;
elseif edgetype=='D'
    structD = getStructD_D;
else
    error('MEMBRANETRACK: invalid edgetype parameter');
end

% compute initial graph from ground truth
[V, E, adjList, faceList] = embryoInitGraph(GT, spacing, structD.optInitSpline, parallel);

% store initial graph structure
VALL = cell(1, N); EALL = cell(1, N);
ADJLIST = cell(1, N); FACELIST = cell(1, N);
VALL{1} = V; EALL{1} = E;
ADJLIST{1} = adjList; FACELIST{1} = faceList;

if verboseG
    mkdir(fname);
    fnamediary = fullfile(fname, sprintf('track_diary.txt'));
    diary(fnamediary);
    f = displayGraph(GT, V, E);
    saveas(f, fullfile(fname, 'initial.png'));
    close(f)
end
% iterate over each image pair
for ii=2:N
    
    disp('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
    fprintf('Tracking images %d -> %d\n', ii-1, ii);
    tic
    
    % compute optimal graph
    I1 = ALL(:,:,ii-1);
    I2 = ALL(:,:,ii);
    
    % get image struct
	structA = struct('I1', {I1}, 'I2', {I2}, 'I2x', {grad(I2)});
    
	% update graph
    structUG = struct('I',{I1},'V',{V},'E',{E},'adjList',{adjList},'faceList',{faceList});
    [ structUG ] = updateGraph(structUG, structC, structD, parallel);
    V = structUG.V;
    E = structUG.E;
    adjList = structUG.adjList;
    faceList = structUG.faceList;
    ADJLIST{ii} = adjList;
    FACELIST{ii} = faceList;
    
    % get graph struct
	structB = struct('V1', {V},'V2',{V},'E1',{E},'E2',{E},'E2new',{E},'adjList',{adjList});
    
	% find optimal graph
	[ V, E, ~ ] = optIterGraph( structA, structB, structC, structD, optOptions);
    VALL{ii} = V;
    EALL{ii} = E;
    
    if verboseG
        f = displayGraph(I2, V, E);
        saveas(f, fullfile(fname, sprintf('result_%d.png', ii)));
        close(f)
    end
	toc
    
end
if verboseG
    diary off
end

data = struct('VALL',VALL,'EALL',EALL,'ADJLIST',ADJLIST,'FACELIST',FACELIST);

end
