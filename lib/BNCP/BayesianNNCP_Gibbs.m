% Bayesian Nonnegative CP Factorization of incomplete tensor
% code by Zecan Yang; Aug. 2022
% % %  All vector is column-wise vectors % % %

classdef BayesianNNCP_Gibbs
    properties
        yTensor;  orgTensor, Map, NoiseT;  % observed tensor, ground_truth tensor; observed tensor index; Noise Tensor;
        rank;  ARD; num_observation; DIMRED;
        alpha_tau_0;  beta_tau_0;
        alpha_lambda_0; beta_lambda_0; 
        dim; order;
        Factors; lambadk; tau; KrFast; nMap; 
        XHat; RMSE_List;
    end
    
    methods
        % create object
        function obj = BayesianNNCP_Gibbs(yTensor, Map, rank, ARD, hyperparameters)
            obj.yTensor = yTensor;
            obj.Map = Map;
            obj.rank = rank;
            obj.ARD = ARD;
            obj.alpha_tau_0 = hyperparameters.alpha_tau_0;
            obj.beta_tau_0 = hyperparameters.beta_tau_0;
            obj.alpha_lambda_0 = hyperparameters.alpha_lambda_0;
            obj.beta_lambda_0 = hyperparameters.beta_lambda_0;
            obj.orgTensor = hyperparameters.orgTensor;
            obj.NoiseT = hyperparameters.NoiseT;
            obj.DIMRED = hyperparameters.DIMRED;
        end

        % initialize parameter
        function self = initialize(self)
            self.dim = size(self.yTensor);        % tensor dimensioan
            self.order = ndims(self.yTensor);     % tensor order

            % request facotr matrix memory
            self.Factors = cell(self.order, 1);   % cell, save Factor Matrix
            for n=1:self.order     
                self.Factors{n} = zeros(self.dim(n), self.rank); 
            end

            self.lambadk = zeros(self.rank, 1);         % init lambad -> ARD
            self.num_observation = sum(self.Map(:));    % observation

            % Initialise lambdak
            for k=1:self.rank
                self.lambadk(k) = self.alpha_lambda_0 / self.beta_lambda_0;
            end

            %  Facotrs Matrix
            for n=1:self.order
                for i=1:self.dim(n)
                    self.Factors{n}(i,:) = Probability_draw.exponential_draw(self.lambadk)';
                end
            end

            % Initialise tau
            self.tau = Probability_draw.gamma_darw(self.alpha_tau(), self.beta_tau());

            initRMSE = lyRMSE(self.orgTensor, double(ktensor(self.Factors)));
            fprintf('Initialise Factors Tensor RMSE: %g\n', initRMSE);
        end


        % % run model
        function self = run(self, RUN_MAX_iterations, AVG_MAX_iterations)
            % for get video
            writerObj = VideoWriter('toy_example_video','Motion JPEG AVI');
            writerObj.FrameRate = 6;  % Set the frames per second
            open(writerObj);

            self.RMSE_List = zeros(RUN_MAX_iterations+1, 1);
            for iter=1:RUN_MAX_iterations
                % update lambdak 
                for k=1:self.rank
                    self.lambadk(k) = Probability_draw.gamma_darw(self.alpha_lambda(), self.beta_lambda(k));
                end

                % update Facotrs Matrix {n}
                for n=1:self.order
                    self = updateCoParameter(self, n);  % KrFast, nMap
                    for r=1:self.rank
                        tau_Factor_n_r = self.tau_Factor(r);
                        mu_Factor_n_r = self.mu_Factor(tau_Factor_n_r, n, r);
                        self.Factors{n}(:,r) = Probability_draw.TN_vector_draw(mu_Factor_n_r, tau_Factor_n_r);
                    end
                end

                % update tau
                self.tau = Probability_draw.gamma_darw(self.alpha_tau(), self.beta_tau());

                % calc error
                self.XHat = double(ktensor(self.Factors));
                self.RMSE_List(iter,1) = lyRMSE(self.orgTensor, double(self.XHat));
                
                if mod(iter, 1)==0
                    fprintf('Iter %3d.  RMSE:%g \n', iter, self.RMSE_List(iter));
                end
                
                %% to plot all information
                if iter==1; gcf = figure('Position',[400 400  300*4 600]); end
                plotALL(self.orgTensor, self.XHat, self.NoiseT, self.yTensor, self.Factors{1}, self.Factors{2}, self.Factors{3}, double(self.RMSE_List(1:iter)));
                %  for save video
                F = getframe(gcf); writeVideo(writerObj, F)

                %% Prune out unnecessary components
                if self.DIMRED==true  && iter >=2
                    Zall = cell2mat(self.Factors);
                    comPower = diag(Zall' * Zall); % comPower'
%                     comTol = sum(self.dim)*eps(norm(Zall,'fro'))*1e-2;
                    comTol = norm(Zall,'fro');
                    rankest = sum(comPower > comTol );
                    if max(rankest)==0
                        disp('Rank becomes 0 !!!');
                        break;
                    end

                    if self.rank ~= max(rankest)
                        indices = comPower > comTol;
                        self.lambadk = self.lambadk(indices);
                        for n=1:self.order
                            self.Factors{n} = self.Factors{n}(:,indices);
                        end
                        self.rank = max(rankest);
                    end
                end
                
            end
            % burn out stage
            self = burn_out(self, AVG_MAX_iterations);
            self.RMSE_List(RUN_MAX_iterations+1,1) = lyRMSE(self.orgTensor, double(self.XHat)); % calc burn-out rmse

            % to plot burn-out information
            plotALL(self.orgTensor, self.XHat, self.NoiseT, self.yTensor, self.Factors{1}, self.Factors{2}, self.Factors{3}, double(self.RMSE_List(1:RUN_MAX_iterations+1)));
            F = getframe(gcf); writeVideo(writerObj, F)
            close(writerObj);
            disp('Successful save video for toy_example!');
        end

        % % burn_out model
        function self = burn_out(self, AVG_MAX_iterations)
            fprintf('AVG—-values： burn-out....................\n');
            Factors_List = cell(self.order, 1); XHat_SUM = zeros(self.dim);
            for n=1:self.order
                Factors_List{n} = zeros(self.dim(n), self.rank); 
            end
            
            for iter=1:AVG_MAX_iterations
                % update lambdak 
                for k=1:self.rank
                    self.lambadk(k) = Probability_draw.gamma_darw(self.alpha_lambda(), self.beta_lambda(k));
                end

                % update Facotrs Matrix {n}
                for n=1:self.order
                    self = updateCoParameter(self, n);  % calc FslashY
                    for r=1:self.rank
                        tau_Factor_n_r = self.tau_Factor(r);
                        mu_Factor_n_r = self.mu_Factor(tau_Factor_n_r, n, r);
                        self.Factors{n}(:,r) = Probability_draw.TN_vector_draw(mu_Factor_n_r, tau_Factor_n_r);
                    end
                    Factors_List{n} = Factors_List{n} + self.Factors{n};
                end

                % update tau
                self.tau = Probability_draw.gamma_darw(self.alpha_tau(), self.beta_tau());

                xhat = double(ktensor(self.Factors));
                XHat_SUM = XHat_SUM + xhat;
               
            end

            for n=1:self.order
                self.Factors{n} = Factors_List{n}./AVG_MAX_iterations;
            end
            self.XHat = XHat_SUM./AVG_MAX_iterations;
            fprintf('Gibbs_burn_out end.....: AVG_XHat_RMSE:%g \n', lyRMSE(self.orgTensor, self.XHat));
        end

        
        %% significant part in our model 
        % % update posterior parameter
        
        % tau_a
        function tau_a = alpha_tau(self)
            tau_a = self.alpha_tau_0 + self.num_observation./2;
        end
        
        % tau_b
        function tau_b = beta_tau(self)
            temp = 0.5*(self.Map .* (self.yTensor-double(ktensor(self.Factors))).^2);
            tau_b = self.beta_tau_0 + sum(temp(:));
        end
        
        % lambda_a
        function lambda_a = alpha_lambda(self)
            lambda_a = self.alpha_lambda_0 + sum(self.dim(:)) ;
        end

        % lambda_b
        function lambda_b = beta_lambda(self,k)
            tmp = 0;
            for n=1:self.order
                tmp = tmp + sum(self.Factors{n}(:,k));
            end
            lambda_b = self.beta_lambda_0 + tmp;
        end

        % Factor Matrix tau index by n
        function res = tau_Factor(self, r)
            tmp1 = self.nMap .* self.KrFast(:,r).^2';
            tmp2 = self.tau .* tmp1;
            res = sum(tmp2, 2);
        end
        
        function res = mu_Factor(self, tau_Facotr, n, r)
            tmp1 = self.yTensor - double(ktensor(self.Factors)) + ktensor_r(self.Factors, r);
            tmp2 = double(tenmat(tmp1, n)) .*  self.KrFast(:, r)';
            tmp3 = self.nMap .* tmp2;
            tmp4 = sum(tmp3, 2); 
            res = 1./(tau_Facotr) .* (-self.lambadk(r) + self.tau .* tmp4);   
        end

        function self = updateCoParameter(self, n)
            self.KrFast = khatrirao_fast(self.Factors{[1:n-1, n+1:self.order]}, 'r');
            self.nMap = double(tenmat(self.Map, n));
        end
    end
end


