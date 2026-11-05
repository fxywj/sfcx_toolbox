classdef nullEngine
    %NULLENGINE Surrogate maps that respect spatial autocorrelation.
    %
    %   kind = 'spin'        spherical rotation, needs parcel coordinates
    %          'msr'         Moran spectral randomisation (eigen sign flips)
    %          'permutation' plain shuffling, anti-conservative

    properties (SetAccess = private)
        kind char
        nPerm double
        seed double
        coh
    end

    properties (Access = private)
        xyz double = []
        lam double = []
        V double = []
    end

    methods
        function obj = nullEngine(coh, kind, nPerm, seed)
            if nargin < 2 || isempty(kind)
                kind = 'msr';
            end
            if nargin < 3 || isempty(nPerm)
                nPerm = 1000;
            end
            if nargin < 4 || isempty(seed)
                seed = 0;
            end
            obj.coh = coh;
            obj.kind = char(kind);
            obj.nPerm = double(nPerm);
            obj.seed = double(seed);

            n = coh.nNodes;
            if strcmp(obj.kind, 'spin')
                if isempty(coh.coords)
                    error('sfcx:null', 'spin null requires coh.coords (nNodes x 3)');
                end
                xyz = double(coh.coords);
                xyz = xyz - mean(xyz, 1);
                r = sqrt(sum(xyz .^ 2, 2));
                r(r < 1e-12) = 1;
                obj.xyz = xyz ./ r;
            elseif strcmp(obj.kind, 'msr')
                a = coh.adjacency();
                d = abs(coh.hierarchy(:) - coh.hierarchy(:)');
                w = a .* exp(-d);
                w(1:n+1:end) = 0;
                m = eye(n) - ones(n) / n;
                wc = m * w * m;
                wc = 0.5 * (wc + wc');
                [obj.V, D] = eig(wc);
                obj.lam = diag(D);
            elseif ~strcmp(obj.kind, 'permutation')
                error('sfcx:null', 'unknown null model ''%s''', obj.kind);
            end
        end

        function S = surrogates(obj, v)
            %SURROGATES nPerm x numel(v) surrogate maps.
            v = double(v(:));
            stream = RandStream('mt19937ar', 'Seed', obj.seed);
            S = zeros(obj.nPerm, numel(v));
            switch obj.kind
                case 'permutation'
                    for i = 1:obj.nPerm
                        S(i, :) = v(randperm(stream, numel(v)))';
                    end
                case 'spin'
                    base = obj.xyz;
                    for i = 1:obj.nPerm
                        R = randomRotation(stream);
                        moved = base * R';
                        d2 = sum(base .^ 2, 2) + sum(moved .^ 2, 2)' - 2 * (base * moved');
                        [~, idx] = min(d2, [], 2);
                        S(i, :) = v(idx)';
                    end
                otherwise
                    proj = obj.V' * v;
                    lam = obj.lam;
                    tol = 1e-6 * (max(abs(lam)) + 1e-12);
                    blocks = {};
                    cur = 1;
                    for k = 2:numel(lam)
                        if abs(abs(lam(k)) - abs(lam(k - 1))) <= tol
                            cur(end+1) = k;
                        else
                            blocks{end+1} = cur;
                            cur = k;
                        end
                    end
                    blocks{end+1} = cur;
                    for i = 1:obj.nPerm
                        signs = ones(size(proj));
                        for b = 1:numel(blocks)
                            blk = blocks{b};
                            if numel(blk) > 1
                                signs(blk) = 2 * (rand(stream, numel(blk), 1) > 0.5) - 1;
                            end
                        end
                        S(i, :) = (obj.V * (proj .* signs))';
                    end
            end
        end
    end
end
