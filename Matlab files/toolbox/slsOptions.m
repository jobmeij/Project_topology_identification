function opt = slsOptions(varargin)
    % Set design parameters of the SLS algorithm
    
    % Tradeoff in speed versus quality, there are 3 options: 'fast',
    % 'best', 'very fast'
    opt.quality = 'fast';
    
    % Maximum number of iterations for efficient estimates
    opt.max_iter = 5;
    
    % When iterations are said to be converged
    opt.T_conv = 1e-1;
    
    % ARX orders to try in Step 1
    opt.ARX_orders = [2 6 10 14 18];
    
    % Weight for Weighted Least Squares    
    opt.WLS = [];
    
    % ARX options, e.g. regularization can be added
    opt.regularization = 'none';
    
    %% Change options according to user specification
    for ii = 1:1:nargin/2
        opt.(varargin{1+2*(ii-1)}) = varargin{2*(ii)};
    end
  
end