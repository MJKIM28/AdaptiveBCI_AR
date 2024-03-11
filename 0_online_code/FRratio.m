function J = FRratio(class1_data,class2_data)
% Combine the data and labels
X = [class1_data; class2_data];
y = [ones(size(class1_data,1), 1); 2*ones(size(class2_data,1), 1)];

% Calculate the overall mean and class-specific means
mu = mean(X);
mu1 = mean(class1_data);
mu2 = mean(class2_data);

% Calculate the within-class scatter matrix
Sw = zeros(size(X, 2));
for i = 1:length(y)
    if y(i) == 1
        Sw = Sw + (X(i,:) - mu1)' * (X(i,:) - mu1);
    else
        Sw = Sw + (X(i,:) - mu2)' * (X(i,:) - mu2);
    end
end
lambda = 1e-6;

Sw = Sw + lambda * eye(size(Sw));
% Calculate the between-class scatter matrix
N1 = sum(y == 1);
N2 = sum(y == 2);
Sb = N1 * (mu1 - mu)' * (mu1 - mu) + N2 * (mu2 - mu)' * (mu2 - mu);

% Calculate the Fisher criterion
% J = det(Sb) / det(Sw);
J = trace(Sb)/trace(Sw);
end