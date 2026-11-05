classdef Cohort
    %COHORT A group of subjects described on a fixed parcellation.
    %
    %   sc, fc are (nSubjects x nEdges) upper-triangle vectors, so every
    %   estimator is a cheap column-wise operation on the same arrays.
    %
    %   coh = sfcx.Cohort(sc, fc, nNodes, networks, hierarchy)
    %   coh = sfcx.Cohort(sc, fc, [], networks, hierarchy, struct('coords', xyz))

    properties
        sc double
        fc double
        bold double = []
        iu double
        iv double
        nNodes double
        nEdges double
        n double
        networks double
        hierarchy double
        coords double = []
        behaviour double = []
        covariates double = []
        state char = 'rest'
        subjectIds double = []
        meta struct = struct()
    end

    methods
        function obj = Cohort(sc, fc, nNodes, networks, hierarchy, opts)
            if nargin < 6 || isempty(opts)
                opts = struct();
            end
            sc = double(sc);
            fc = double(fc);
            if isvector(sc)
                sc = sc(:)';
                fc = fc(:)';
            end
            if ~isequal(size(sc), size(fc))
                error('sfcx:Cohort', 'sc and fc must have the same size, got [%s] and [%s]', ...
                    num2str(size(sc)), num2str(size(fc)));
            end
            if ndims(sc) ~= 2
                error('sfcx:Cohort', 'sc and fc must be 2-D (subjects x edges)');
            end
            [obj.n, obj.nEdges] = size(sc);
            if isempty(nNodes)
                nNodes = round((1 + sqrt(1 + 8 * obj.nEdges)) / 2);
            end
            obj.nNodes = double(nNodes);
            if obj.nEdges ~= obj.nNodes * (obj.nNodes - 1) / 2
                error('sfcx:Cohort', 'nEdges=%d is not the upper triangle of nNodes=%d', ...
                    obj.nEdges, obj.nNodes);
            end
            obj.sc = sc;
            obj.fc = fc;
            [obj.iu, obj.iv] = upperTriIndex(obj.nNodes);
            if nargin < 4 || isempty(networks)
                networks = zeros(1, obj.nNodes);
            end
            if nargin < 5 || isempty(hierarchy)
                hierarchy = linspace(1, 7, obj.nNodes);
            end
            obj.networks = double(networks(:)');
            obj.hierarchy = double(hierarchy(:)');
            if numel(obj.networks) ~= obj.nNodes
                error('sfcx:Cohort', 'networks must have nNodes entries');
            end
            if numel(obj.hierarchy) ~= obj.nNodes
                error('sfcx:Cohort', 'hierarchy must have nNodes entries');
            end
            if isfield(opts, 'coords') && ~isempty(opts.coords)
                if ~isequal(size(opts.coords), [obj.nNodes, 3])
                    error('sfcx:Cohort', 'coords must be (nNodes, 3)');
                end
                obj.coords = double(opts.coords);
            end
            if isfield(opts, 'behaviour') && ~isempty(opts.behaviour)
                b = double(opts.behaviour(:));
                if numel(b) ~= obj.n
                    error('sfcx:Cohort', 'behaviour must have nSubjects entries');
                end
                obj.behaviour = b;
            end
            if isfield(opts, 'covariates') && ~isempty(opts.covariates)
                c = double(opts.covariates);
                if isvector(c)
                    c = c(:);
                end
                if size(c, 1) ~= obj.n
                    error('sfcx:Cohort', 'covariates must have nSubjects rows');
                end
                obj.covariates = c;
            end
            if isfield(opts, 'state') && ~isempty(opts.state)
                obj.state = char(opts.state);
            end
            if isfield(opts, 'subjectIds') && ~isempty(opts.subjectIds)
                obj.subjectIds = double(opts.subjectIds(:));
            else
                obj.subjectIds = (1:obj.n)';
            end
            if isfield(opts, 'bold') && ~isempty(opts.bold)
                b = double(opts.bold);
                if ismatrix(b)
                    b = reshape(b, 1, size(b, 1), size(b, 2));
                end
                if ndims(b) ~= 3
                    error('sfcx:Cohort', 'bold must be 3-D (nSubjects, nNodes, nTime)');
                end
                if size(b, 1) ~= obj.n || size(b, 2) ~= obj.nNodes
                    error('sfcx:Cohort', 'bold must be (%d, %d, nTime), got [%s]', ...
                        obj.n, obj.nNodes, num2str(size(b)));
                end
                obj.bold = b;
            end
        end

        function k = nNetworks(obj)
            k = max(obj.networks) + 1;
        end

        function a = adjacency(obj, threshold)
            %ADJACENCY Group-mean structural matrix, binarised.
            if nargin < 2 || isempty(threshold)
                threshold = 0;
            end
            m = obj.groupMatrix('sc');
            a = double(m > max(threshold, 1e-9));
        end

        function m = groupMatrix(obj, which)
            %GROUPMATRIX Group-mean connectivity as a symmetric node matrix.
            if nargin < 2
                which = 'sc';
            end
            if strcmpi(which, 'fc')
                v = nanMeanVec(obj.fc, 1);
            else
                v = nanMeanVec(obj.sc, 1);
            end
            m = vecToMat(v, obj.iu, obj.iv, obj.nNodes);
        end

        function c = regionalCoupling(obj)
            %REGIONALCOUPLING c_i = corr(SC_i., FC_i.) per subject and region.
            c = nan(obj.n, obj.nNodes);
            for i = 1:obj.nNodes
                idx = find((obj.iu == i) | (obj.iv == i));
                if numel(idx) < 3
                    continue;
                end
                a = double(obj.sc(:, idx));
                b = double(obj.fc(:, idx));
                a = a - mean(a, 2);
                b = b - mean(b, 2);
                sa = sqrt(sum(a .^ 2, 2));
                sb = sqrt(sum(b .^ 2, 2));
                ok = (sa > 1e-12) & (sb > 1e-12);
                val = nan(obj.n, 1);
                val(ok) = sum(a(ok, :) .* b(ok, :), 2) ./ (sa(ok) .* sb(ok));
                c(:, i) = val;
            end
        end

        function [a, b] = splitHalf(obj, seed)
            %SPLITHALF Two disjoint index sets of equal size.
            if nargin < 2
                seed = 0;
            end
            rng = RandStream('mt19937ar', 'Seed', seed);
            perm = randperm(rng, obj.n);
            h = floor(obj.n / 2);
            a = perm(1:h);
            b = perm(h+1:2*h);
            a = a(:);
            b = b(:);
        end
    end
end
