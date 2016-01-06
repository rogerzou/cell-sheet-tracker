function [ V2, E2, c ] = optIterGraph( structA, structB, structC, structD, optOptions )
%OPTITERGRAPH Compute the optimal graph between two images.
% Given two images I1, I2, where V1 is the vertex set for I1, V2 is the
% corresponding vertex set for I2, and E1 are the fixed splines in I1
% and E2 are corresponding splines in I2, find the optimal placement
% of points in V2 with correspondigng i2, that minimizes the error
% between the graph cost on I1 and I2. Optimization performed using Newton's
% Method.
% 
% INPUTS
% structA: images struct.
%   - I1: (R^2 matrix) image I1, the previous image in a sequence.
%   - I2: (R^2 matrix) image I2, the current image in a sequence.
%   - I2x: (2 x 1 struct of R^2 matrix) gradient of image I2.
% structB: graph struct.
%   - V1: (cell of (d x 1) vectors), vertex set in image I1.
%   - V2: (cell of (d x 1) vectors), initial vertex set in image I2.
%   - E1: (cell of spline structs), edge set in image I1.
%   - E2: (cell of spline structs), initial edge set in image I2.
%   - E2new: (cell of structs) the splines in image I2, some of which are
%   adjacent to vertex v in image I2 and moved in the CURRENT ITERATION.
%   If E2new is set equal to E2, then that assumes all vertices
%   are updated simultaneously in Newton's method.
%   - adjList
% structC: parameters struct.
%   - alpha: weight of vertex cost (1 = all vertex, 0 = all edge)
%   - l: (odd scalar) width of the spline edge in # of pixels.
%   - w: (odd scalar) width of the vertex window in # of pixels.
% structD: functions struct.
%   - fvert
%   - fedge
% optOptions: options struct with fields: 
%   - parallel
%   - verboseE
%   - verboseG
%
% OUTPUTS
% V2: (cell of (d x 1) vectors), optimal vertex set in image I2.
% E2: (cell of spline structs), optimal edge set in image I2.
% c: (scalar) cost of optimal graph.
%
% CALLEE functions
%   optIterSpline
%   Newton
%   integrate
%   D2gcostDv2
%   jacobian
%   vNeighborE
%   splineEvalEven
%
% @author Roger Zou
% @date 8/16/15

% get and validate input variables
V20 = structB.V2;
E1 = structB.E1;
E20 = structB.E2;
adjList = structB.adjList;
l = structC.l;
w = structC.w;
updateSplineEndPts = structD.updateSplineEndPts;
d = 2;
if ~mod(l, 2)
    error('OPTITERGRAPH: l must be odd');
end
if ~mod(w, 2)
    error('OPTITERGRAPH: w must be odd');
end

% parse optimization parameters
parallel = false;
verboseE = false;
verboseG = false;
siftflow = false;
if nargin==nargin('optIterGraph') && ~isempty(optOptions)
    if any(strcmp('parallel',fieldnames(optOptions)))
        parallel = optOptions.parallel;
    end
    if any(strcmp('verboseE',fieldnames(optOptions)))
        verboseE = optOptions.verboseE;
    end
    if any(strcmp('verboseG',fieldnames(optOptions)))
        verboseG = optOptions.verboseG;
    end
    if any(strcmp('siftflow',fieldnames(optOptions)))
        siftflow = optOptions.siftflow;
    end
end

if verboseG
    fprintf('\n############ OPTITERGRAPH START ############\n');
end

% pack initial vertex set V0, edge set SplinesJ0 into a struct
X0 = pack(V20, E20);

% SIFT flow to compute initial vertex positions
if siftflow
    step = SIFTflow(structA.I1, structA.I2, V20);
    X0 = graphStep(X0, step);
end

% Compute optimal graph with Newton's method
[X, c] = Newton(@graphCost, @graphOptA, @graphStep, X0, verboseG);

% unpack X to retrieve final vertex and edge set
[V2, E2] = unpack(X);

if verboseG
    fprintf('\n############ OPTITERGRAPH END ############\n');
end


    %% HELPER FUNCTIONS

    function [neVI] = ismergedV(V)
        neVItmp = ~sum(isnan(cell2mat(V')))';
        neVI = false(numel(neVItmp)*d, 1);
        neVI(1:d:end-1) = neVItmp;
        neVI(2:d:end) = neVItmp;
    end
    
    % obtain struct X of vertex and edge set for input to optimization
    function X = pack(V, E)
        X = struct;
        X.vertices = V;
        X.edges = E;
    end

    % obtain vertex set and splines from struct X
    function [V, E, N, M, n] = unpack(X)
        V = X.vertices;
        E = X.edges;
        if nargout > 2
            N = length(V);
            M = length(E);
            n = N*d;    % number of optimization variables
        end
    end

    % compute the graph cost c, with $0 \le \alpha \le 1$ as the weight on
    % the vertex cost. Then the weight on edge cost is (1-alpha). n is the
    % number of optimization variables
    function [c, n] = graphCost(X)
        [V2k, E2k, ~, ~, n] = unpack(X);
        structBk = structB;
        structBk.V2 = V2k;
        structBk.E2 = E2k;
        structBk.E2new = E2k;
        c = gcost(structA, structBk, structC, structD, parallel);
    end

    %%% Compute analytic gradient and Gramian for Taylor approximation %%%
    function [g, Gr] = graphOptA(X)
        [V2k, E2k] = unpack(X);
        structBk = structB;
        structBk.V2 = V2k;
        structBk.E2 = E2k;
        structBk.E2new = E2k;
        [g, Gr] = D2gcostDv2(structA, structBk, structC, structD, parallel);
        neVI = ismergedV(V2k);
        g = g(neVI);
        Gr = Gr(neVI,neVI);
    end

    % newton step - move vertices, find optimal new splines
    function X = graphStep(X, step)
        
        [V2k, E2k, N, M] = unpack(X);
        neVI = ismergedV(V2k);
        steptmp = NaN(N*d,1);
        steptmp(neVI) = step;
        step = steptmp;
        
        % store all endpoint info in updateSI
        updateSIa = NaN(d,M);
        updateSIb = NaN(d,M);
        for ii=1:M
            E2ki = E2k{ii};
            if ~isempty(E2ki)
                sctrl = E2ki.control;
                updateSIa(:,ii) = sctrl(:,1);
                updateSIb(:,ii) = sctrl(:,end);
            end
        end
        
        % iterate over each vertex
        r1 = 1 + d * ( (1:N) - 1 );
        r2 = d * (1:N);
        for ii=1:N
            
            % get initial vertex position
            vk0 = V2k{ii};
            
            if ~sum(isnan(vk0))
                % get newton step for this vertex
                vstep = step(r1(ii):r2(ii));
                % update vertex set
                V2k{ii} = vk0 + vstep;

                % get indices of edges adjacent to initial vertex position
                [eindices] = vNeighborE(ii, V2k, E2k, adjList);

                % iterate over the splines attached to the shifted vertex
                for jj=eindices'

                    % get the current spline endpoints
                    a = updateSIa(:,jj);
                    b = updateSIb(:,jj);

                    % take the newton step for the vertex of interest
                    if a(1)==vk0(1) && a(2)==vk0(2)
                        a = a + vstep;
                    elseif b(1)==vk0(1) && b(2)==vk0(2)
                        b = b + vstep;
                    end

                    % update the spline endpoints
                    updateSIa(:,jj) = a;
                    updateSIb(:,jj) = b;

                end
            end
        end

        % get all optimal splines with new endpoints (PARFOR POSSIBLE)
        if parallel
            parfor ii=1:M
                
                E20i = E20{ii};
                if ~isempty(E20i)
                    % update the spline with the new endpoints
                    sV = updateSplineEndPts( E20i, updateSIa(:,ii), updateSIb(:,ii) );

                    % recompute optimal spline
                    [sV, ~] = optIterSpline( E1{ii}, sV, ...
                        structA, structC, structD, verboseE);
                    E2k{ii} = sV;
                end
                
            end
        else
            for ii=1:M
                
                E20i = E20{ii};
                if ~isempty(E20i)
                    % update the spline with the new endpoints
                    sV = updateSplineEndPts( E20i, updateSIa(:,ii), updateSIb(:,ii) );

                    % recompute optimal spline
                    [sV, ~] = optIterSpline( E1{ii}, sV, ...
                        structA, structC, structD, verboseE);
                    E2k{ii} = sV;
                end
                
            end
        end
        
        X = pack(V2k, E2k);
        
    end

end
