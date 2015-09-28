function [mu,Sigma] = updateStatistics_Lacki15(mu,Sigma,theta,i,X,XX,r)                       

gamma = 1/i;
mu = (1-gamma)*mu + gamma*theta;
Sigma = (1-gamma)*Sigma + gamma*(theta-mu)*(theta-mu)';

% Regularisation
[~,p] = cholcov(Sigma,0);
if p ~= 0
    Sigma = Sigma + r*eye(size(dtheta));
end

end



