function combine(string,n)
% filename = input('Enter the string:  ', 's');

comb=zeros(n,41); c=zeros(2,41);

for i=1:n
    load ([string num2str(i)]);
    corr_var1_var2(i,:) = c(2,:);
%     clear c
end
x_lag = c(1,:);

y_data = reshape(corr_var1_var2', 1, prod(size(corr_var1_var2)));
x_data = repmat(x_lag, 1, size(corr_var1_var2,1));
[max_x_residuals] = spline_fit_bootstrap(x_data, y_data);

end