% inputs
budgetSol.x = 1e6 *[4.6360 0.6364 0.5035 4.2241]';
boxSol.x = 1e6 * [1.6193 0.6543 0.5547 7.1716]';
nomSol.x = 1e6 * [0.2345 0.0359 0.0167 9.7129]';
n = 4;
iter = 10000;
abar = [0.608 0.386 0.430 0.864]';
ahat = [0.257 0.114 0.127 0.329]';
cs = [8552.695 7195.047 6045.412 4381.518]';
ds = [10000 5314.612 10000 10000]';
sum_h = zeros(n,1)
h_save = zeros(n,iter)
%get the validation of samples
for i = 1:iter
    %z = get_z(n);
    z = rand(n, 1);
    a = abar.*(1 - ahat .* z);
    h_save(:,i) = cs.*(1+boxSol.x./ds).^a - cs;
    %h = cs.*(1+boxSol.x./ds).^a - cs;
    %h_save(i) = cs.*(1+nomSol.x./ds).^a - cs;
end

sum_h = sum(h_save);
hist(sum_h);
%boxplot(h_save');

% define the uncertainty
function z = get_z(n)
x = rand(n, 1); % Start with 3 random numbers that don't sum to 1.
z = x / sum(x)  % Normalize so the sum is 1.
end


