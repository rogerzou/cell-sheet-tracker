function [ structUG ] = updateGraph( structUG, structC, structD, parallel)
%UPDATEGRAPH Update the graph topology between each optimization run.
% - Merge vertices that are too close and merge associated edges.
% - Update the number of internal control points.
%
% INPUTS
% structUG: updategraph struct.
%   - I: (2D matrix) image associated with V and E.
%   - V: (cell of (d x 1) vectors) vertex set.
%   - E: (cell of spline structs) edge set.
%   - adjList: (Nx1 cell array): for index i, contains a (2xK) matrix of K
%   vertex-edge index pairs of vertices and edges adjacent to vertex i.
%   - faceList: (Fx1 cell array) for index i, contains a (2xL) matrix of L
%   vertex,edge indices of vertices and edges that are part of cell i.
% structC: parameters struct.
%   - interval: (scalar) minimum distance between vertices before merging.
%   - order: (integer) the order of each spline.
%   - spacing: (scalar) spacing between each control point in pixels.
% structD: functions struct.
% parallel: (boolean) true to activate parallelized code.
%
% OUTPUTS
% structUG: updated updategraph struct.
%
% CALLEE functions
%   uvEdge
%   vNeighborE
%   updateSplineEndPoints
%
% @author Roger Zou
% @date 8/15/15

% validate inputs
if nargin==nargin('updateGraph')-1;
    parallel = false;
end

% get input variables
I = structUG.I;
V = structUG.V;
E = structUG.E;
adjList = structUG.adjList;
faceList = structUG.faceList;
interval = structC.interval;
spacing = structC.spacing;
optInitSpline = structD.optInitSpline;
updateSplineEndPts = structD.updateSplineEndPts;
order = 3;
d = 2;

%% Merge vertices that are too close and merge associated edges.

% save copy of original input splines
Eold = E;

counter = 0;
while true

    % get sizes
    N = numel(V);

    % sort the vertices by distance
    Varray = NaN(N, d);
    for ii=1:N
        Varray(ii,:) = V{ii}';
    end
    D = triu( squareform(pdist(Varray)), 1);
    D(tril(true(N))) = Inf;
    [Vals, Ind] = sort(D(:));
    Vals(Vals>interval) = [];
    Vals(isnan(Vals)) = [];
    
    % if vertex set not empty and distance greater than interval
    if isempty(Vals)
        break
    else
        
        % get the two vertices with smallest pairwise distance
        [U1, U2] = meshgrid(1:N, 1:N);
        U1 = U1(:); U2 = U2(:);
        uaInd = U1(Ind(1));
        ubInd = U2(Ind(1));
        ua = V{uaInd};
        ub = V{ubInd};

        % get indices of edges adjacent to both ua and ub
        [edgeIa] = vNeighborE(uaInd, V, E, adjList);
        [edgeIb] = vNeighborE(ubInd, V, E, adjList);

        % get indices of edges e between ua and ub
        eInds = intersect(edgeIa, edgeIb);

        % update indices of edges adjacent to ua, ub, minus e
        edgeIa = setdiff(edgeIa, eInds);
        edgeIb = setdiff(edgeIb, eInds);

        % make u_new the center of ua and ub
        unew = (ua + ub) ./ 2;

        % make edges formerly adjacent to ua, to unew
        if ~isempty(edgeIa)
            for edgei=edgeIa'
                s = E{edgei};
                sa = s.control(:,1);
                sb = s.control(:,end);
                if ua==sa
                    va = unew;
                    vb = sb;
                elseif ua==sb
                    va = sa;
                    vb = unew;
                else
                    error('UPDATEGRAPH: unexplained error!');
                end
                E{edgei} = updateSplineEndPts(s, va, vb);
            end
        end

        % make edges formerly adjacent to ub, to unew
        if ~isempty(edgeIb)
            for edgei=edgeIb'
                s = E{edgei};
                sa = s.control(:,1);
                sb = s.control(:,end);
                if ub==sa
                    va = unew;
                    vb = sb;
                elseif ub==sb
                    va = sa;
                    vb = unew;
                else
                    error('UPDATEVERTICES: unexplained error!');
                end
                E{edgei} = updateSplineEndPts(s, va, vb);
            end
        end

        % move ua to unew
        V{uaInd} = unew;
        
        % delete ub
        V{ubInd} = [NaN; NaN];
        
        % delete e
        if ~isempty(eInds)
            for ii=1:numel(eInds)
                E{eInds(ii)} = [];
                Eold{eInds(ii)} = [];
            end
        end
        
        % increment counter
        counter = counter + 1;
        
        % update adjList
        adjListOld = adjList;
        for ii=1:numel(adjListOld);
            adjListi = adjListOld{ii};
            if ii==uaInd
                adjListi(:, adjListi(1,:)==ubInd) = [];
                adjListb = adjListOld{ubInd};
                adjListb(:, adjListb(1,:)==uaInd) = [];
                adjListi = [adjListi, adjListb]; %#ok<AGROW>
            else
                adjListi(1, adjListi(1,:)==ubInd) = uaInd;
            end
            adjList{ii} = adjListi;
        end
        adjList{ubInd} = zeros(2,0);
        
        % update faceList
        for ii=1:numel(faceList)
            faceListi = faceList{ii};
            if ~isempty(faceListi)
                faceListv = faceListi(1,:);
                faceListe = faceListi(2,:);
                if ismember(eInds, faceListe)   % if edges are part of cell
                    faceListv(faceListv==ubInd) = 0;    % set deleted face vertex components to zero (bc UINT32 TYPE)
                    for jj=1:numel(eInds)
                        faceListe(faceListe==eInds(jj)) = 0;    % set deleted face edge components to zero (bc UINT32 TYPE)
                    end
                else                            % if edges are not part of cell
                    faceListv(faceListv==ubInd) = uaInd;
                end
                if numel(faceListv)==1  % delete face if only one vertex
                    faceList{ii} = [];
                else                    % put current face back into structure
                    faceListi = [faceListv; faceListe];
                    faceList{ii} = faceListi;
                end
            end
        end
    end
end

% display info
% fprintf('\nUPDATEGRAPH: %d vertices merged.\n', counter);

% get nonempty splines
enempty = ~cellfun('isempty', E);
neE = E(enempty);
neEold = Eold(enempty);

% find optimal splines (PARFOR POSSIBLE)
M = numel(neE);
structA = struct('I1', I, 'I2', I, 'I2x', {grad(I)});
if parallel
    parfor ii=1:M
        neE{ii} = optIterSpline(neEold{ii}, neE{ii}, ...
            structA, structC, structD, false);
    end
else
	for ii=1:M
        neE{ii} = optIterSpline(neEold{ii}, neE{ii}, ...
            structA, structC, structD, false);
	end
end
E(enempty) = neE;


%% Update the number of internal control points.

% get nonempty splines
enempty = ~cellfun('isempty', E);
neE = E(enempty);

% get number of splines
M = numel(neE);

% recompute the control points for each spline (PARFOR POSSIBLE)
if parallel
    parfor ii=1:M
        scurve = neE{ii}.curve;
        % dynamically compute number of control points
        k = nCtrlPts(scurve, spacing);
        % create optimal spline from sampling
        neE{ii} = optInitSpline(scurve, order, k);
    end
else
    for ii=1:M
        scurve = neE{ii}.curve;
        % dynamically compute number of control points
        k = nCtrlPts(scurve, spacing);
        % create optimal spline from sampling
        neE{ii} = optInitSpline(scurve, order, k);
    end
end
E(enempty) = neE;

%% repackage for output
structUG.E = E;
structUG.V = V;
structUG.adjList = adjList;
structUG.faceList = faceList;

