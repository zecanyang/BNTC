% %[Author Notes]%   
%   Author: Zecan Yang
%   Email : zecanyang@gmail.com
%   Aug. 2022

classdef Probability_draw
    methods(Static)

        % Gamma function
        function res = gamma_darw(alpha, beta)
            shape = alpha;
            scale = 1./beta;
            res = gamrnd(shape, scale);
        end

        % Exponential function
        function res = exponential_draw(hyperparam)
            scale = 1.0 ./ hyperparam;
            res = exprnd(scale);
        end

        % Truncated Gaussian function
        function draws = TN_vector_draw(mus, taus)
            sigmas = 1./ sqrt(taus);
            draws = zeros(length(mus), 1);

            for i=1:length(mus)
                if taus(i) == 0 
                    draws(i) = 0; 
                else
                    tmp = rtnorm(0, inf, mus(i), sigmas(i)); %  Truncated singular value
                    if tmp<0 || tmp==inf || isnan(tmp)       %  Recheck
                        draws(i) = 0;
                    else
                        draws(i) = tmp;
                    end
                end
            end
        end

        function res = TN_vector_expection(mus,taus)
            sigmas = 1./ sqrt(taus);
            res = zeros(length(mus), 1);
            x = -mus./sigmas;
            lambdax = normpdf(x)./(0.5*erfc(x./sqrt(2)));
            expectation = mus + sigmas .* lambdax;
            for i=1:length(expectation)
                if mus(i) < -30 * sigmas(i)
                    expectation(i) = 1./(abs(mus(i))*taus(i));
                else
                    continue;
                end
            end
            for i=1:length(expectation)
                if(expectation(i)>=0.0 && expectation(i)~=inf && expectation(i)~=-inf && ~(isnan(expectation(i))))
                    res(i) = expectation(i);
                else
                    res(i) = 0.;
                end
            end
        end
        
        function res = TN_vector_variance(mus,taus)
            sigmas = 1./ sqrt(taus);
            res = zeros(length(mus), 1);
            x = -mus./sigmas;
            lambdax = normpdf(x)./(0.5*erfc(x/sqrt(2)));
            deltax = lambdax.*(lambdax-x);
            variance = sigmas.^2 .* (1 - deltax);
            for i=1:length(variance)
                if mus(i)< -30*sigmas(i)
                    variance(i) = (1./(abs(mus(i))*taus(i))).^2;
                else
                    continue;
                end
            end
            for i=1:length(variance)
                if(variance(i)>=0.0 && variance(i)~=inf && variance(i)~=-inf && ~(isnan(variance(i))))
                    res(i) = variance(i);
                else
                    res(i) = 0.;
                end
            end
        end
    end
end