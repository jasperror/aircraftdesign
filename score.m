% ENGR 493 Intro to Aircraft Design
% Project Phase 2: Wing Design - Scoring
% Jasper Palmer
classdef score
    properties (Constant)
        AoAstall = 0.2;
        Ld = 0.8;
    end
    methods (Static)
        % Trade studies:
        % - Ld vs alpha_stall
        % - 
        % best b, lambda, Lambda_LE combination
        % goal: highest LD_cr with reasonable tradeoffs
        
        function score = score(Ld, AoA_stall)
            score = Ld.*(score.Ld/max(Ld,[],'all')) + AoA_stall.*(score.AoAstall/max(AoA_stall,[],'all'));
        end
    end
end