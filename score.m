classdef score
    % Calculates the score of each alternative in order to allow for design
    % optimization
    
    properties % Inputs
        stealth % Somehow calculate the stealth score (or split it up)
        cost % Aircraft cost, CA$
        perf % Aircraft performance
        % M_dash_lo % Strike dash mach, []
        % M_dash_cr % Dogfight dash mach, []
        % psidot_sturn % Dogfight sustained turn rate, deg/s
        % range % Range, m
        % t_lo % Dogfight loiter time, s
        % N_z % Maximum vertical load factor, []
        % m_d % Drop weight, kg
    end

    properties (Constant)
        cm % Comparison matrices
        norm = [1, 3, 1]; % Normalization methods (1 for benefit, 3 for cost, 5 for nonmonotonic)
        nonm = [0, 0, 0]; % Nonmonotomic ideal values
    end

    properties (Dependent) % Outputs
        z % Scores
        dm % Decision matrix
        ndm % Normalized decision matrix
        w % Criteria weights
    end
    
    methods
        function z = get.z(obj)
            % Calculate both the linear and vector normalization scores
            z.lin = 1;
            z.vec = 1;
        end
        function dm = get.dm(obj)

        end
        function ndm = get.ndm(obj)
            % Create the normalized decision matrix
            for j = 1:size(obj.dm,2)
                switch normdict(j)
                    case 1
                        ndm(:,j) = linear_ben(obj.dm(:,j));
                    case 2
                        ndm(:,j) = vec_ben(obj.dm(:,j));
                    case 3
                        ndm(:,j) = lin_cost(obj.dm(:,j));
                    case 4
                        ndm(:,j) = vec_cost(obj.dm(:,j));
                    case 5
                        ndm(:,j) = nonmonotonic(obj.dm(:,j), obj.nonm(j));
                end
            end
        end
        function w = get.w(obj) % Split up into a main obj.w which depends on a constant which obj.w_comparisonmatrix or obj.w_entropy, etc.?
            % Calculate weights from comparison matrices (or whatever we
            % end up going with)
            w = 1;
        end
    end

    methods (Static) % Normalization Methods
        function c = linear_ben(col)
            c = col./max(col);
        end
        function c = vec_ben(col)
            c = col./sqrt(sum(col.^2));
        end
        function c = lin_cost(col)
            c = (1./col)./max(1./col);
        end
        function c = vec_cost(col)
            c = (1./col)./sqrt(sum((1./col).^2));
        end
        function c = nonmonotonic(col, ideal)
            z = (col-ideal)./std(col);
            c = exp(-0.5.*z.^2);
        end
    end
end

