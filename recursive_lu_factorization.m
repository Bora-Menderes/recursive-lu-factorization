%--------Sherman's March Factorization for A = LU-------------

clc

%-----FUNDAMENTAL FUNCTIONS------

%function to apply forward_substitution, Lx = b
function solution = forward_substitution(L, b) %L lower triangular, b vector
    dim = size(L, 1); %get the dimension of L

    x = zeros(dim, 1); %initialize a vector x
    x(1) = b(1) / L(1, 1);
    for i = 2:dim
        x(i) = (b(i) - L(i, 1:i-1) * x(1:i-1)) / L(i, i); %forward subsitution
    end
    solution = x;
end

%function to apply backward substitution, Ux = b
function solution = backward_substitution(U, b)
    dim = size(U, 1); %get the size of U

    x = zeros(dim, 1); %initialize a vector x for the solution
    x(dim) = b(dim) / U(dim, dim);
    for i = dim-1:-1:1
        x(i) = (b(i) - U(i, i+1:dim) * x(i + 1:dim)) / U(i, i); %backwards subsitution
    end
    solution = x;
end

%function to apply LU Factorization with Sherman's March
function [L, U] = factorization(A)
    dim = size(A, 1); % get the dimension of A

    %initialize L and U with base cases
    L = zeros(dim, dim);
    U = zeros(dim, dim);
    L(1, 1) = 1; U(1, 1) = A(1, 1);

    for i = 2:dim
        L(i, i) = 1; %set the diagonal entry of L to 1

        %Compute new row vector for L A = LU
        L(i, 1:i-1) = forward_substitution(U(1:i-1, 1:i-1)', A(i, 1:i-1)')'; 
        %Note: Must use transposes to align vectors

        %Compute the new column vector for U
        U(1:i-1, i) = forward_substitution(L(1:i-1, 1:i-1), A(1:i-1, i));

        %Compute the new diagonal entry for U
        U(i, i) = A(i, i) - L(i, 1:i-1) * U(1:i-1, i);
    end
end

%----Time Recording------
%initialize an array to store timestamps
%rows signify the dimension, columns represent different tries with same
%dimension
recordedtimes = zeros(200, 30); 

%Record in time taken to calculate LU
for n = 1:200
    A = hilb(n);
    for i = 1:30
        tic
        [L, U] = factorization(A);
        toc
        recordedtimes(n, i) = toc;
    end
end

%set up the confidence interval

meanTimes = mean(recordedtimes, 2); %take the averages
stdTimes  = std(recordedtimes, 0, 2); %compute standard dev.

%calculate the %95 confidence interval
Nsamples = 30; 
CI95 = 1.96 * stdTimes / sqrt(Nsamples);

upper = meanTimes + CI95;
lower = meanTimes - CI95;

%plot 

n = 1:200;

figure; hold on;

%Draw the mean curve
plot(n, meanTimes, 'b', 'LineWidth', 2);

%Apply the confidence interval, use shades
fill([n fliplr(n)], [upper' fliplr(lower')], ...
     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

%set up plot labels
xlabel('Matrix dimension n');
ylabel('Runtime (in seconds)');
title("Runtime of Sherman's March with 95% Confidence Intervals");
legend('Mean runtime', '95% confidence band');


%---Relative Reproduction Error----

%create a 200x2 matrix to store errors for each n = 1,2,...,200
%1st column represents errors from our Sherman's March factorization
%2nd column represents errors from MATLAB's lu() factorization
relative_production_errors = zeros(200, 2);

for n = 1:200
    A = hilb(n);

    [L1, U1] = factorization(A); %factorize A with Sherman's
    [L2, U2, P] = lu(A); %factorize A with matlab lu()

    %Relative Error For Sherman's
    abs_error_sherman = norm(A - L1*U1, "fro"); %First, calculate the absolute error
    %calculate and store the relative error
    relative_production_errors(n, 1) = abs_error_sherman / norm(A, "fro");

    %Relative Error For MATLAB's lu()
    abs_error_matlab_lu = norm(P*A - L2*U2, "fro"); %First, calculate the absolute error
    %calculate and assign the relative production error
    relative_production_errors(n, 2) = abs_error_matlab_lu / norm(A, "fro");
end

%-----FINANCIAL APPLICATION-------

%K - 100 = 150

function I = minimum_variance_portfolio(C,M)

    N = size(C, 1); %determine the dimension

    vector_one = ones(N, 1); %Create a column vectors of 1
    I = zeros(N, 2); %Matrix to store I_sherman and I_lu in its
    %columns
   
    %Note: C is 20x250, rows are assets
    %columns are days; C' is the other way
    cov_matrix = cov(C');
    %Factorie Cov(C') into LU to solve appropriate systems
    [L1, U1] = factorization(cov_matrix); %Sherman's
    [L2, U2, P] = lu(cov_matrix); %MATLAB
   
    %For Sherman's Charge: calculate I_sherman
    %Solve covariance*x = L1*U1*x = 1 // L1*y = 1, U1x = y
    y = forward_substitution(L1, vector_one);
    x1 = backward_substitution(U1, y);

    w1 = x1 / (vector_one' * x1);
    for i=1:N
        I(i, 1) = w1(i) * M;
    end

    %For MATLAB's lu(): calculate I_lu
    %Solve covariance*x = 1, i.e., L2*U2*x = P * 1 // L2*y=P * 1, U2*x = y
    y = forward_substitution(L2, P * vector_one);
    x2 = backward_substitution(U2, y);

    w2 = x2 / (vector_one' * x2); 

    for i=1:N %store invesment of each asset
        I(i, 2) = w2(i) * M;
    end
end

function profit = calculate_profit(assetamount, C, M)

    profit = zeros(250, 2); %store profits for each j > K - 100. To plot more
    %smoothly this profit matrix has 250 rows, but will only use to last
    %100
    invesment_of_assets = minimum_variance_portfolio(C(1:assetamount,1:150), M);
    
    %find prices of assets at day 150, and buy these assets. Do this by
    %updating (overwriting) previous invesment_of_assets(i) so that it becomes the number
    %of stocks that are bought. 

    for i=1:assetamount
        for j=1:2
            invesment_of_assets(i, j) = invesment_of_assets(i, j) / C(i, 150);
        end
    end

    %Finally calculate profit for each day j > K - 100
    for j = 150:250
            for i=1:2
                profit(j, i) = invesment_of_assets(1:assetamount,i)' * C(1:assetamount, j) - M;
            end
    end
end






