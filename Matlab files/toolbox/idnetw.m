function net = idnetw(varargin)
    %IDNETW Constructs a polynomial network model
    %
    % Use it with pre-programmed network structures:
    % net = idnetw(structure,L,K,order)
    % 
    % Add the topology of G:
    % net = idnetw(structure,L,K,order,topoG)
    % 
    % With manually specified structures:
    % net = idnetw(topo,n)
    %
    % Or add fixed dynamics:
    % net = idnetw(topo,n,fixed)
    %
    % structure is a string: 
    %   'fixedR' models R=I of dimension L by K, there are no parameters in R. 
    %            The H is diagonal, and does contain parameters;
    %   'fixedRfullH' models R=I of dimension L by K, there are no parameters in R. 
    %                 The H is full, i.e. all elements are non-zero and contain parameters.		
    % 	'diagR' models R as a diagonal matrix of dimension L by K, and it does contain parameters.
    %           The H is diagonal, and does contain parameters;
    % 	'diagRfullH' models R as a diagonal matrix of dimension L by K, and it does contain parameters.
    %                The H is full, i.e. all elements are non-zero and contain parameters.
    %
    % If you want to add the structure of $G$ to one of these pre-specified model structures, then the additional argument \texttt{adjG} can be added. 
    % \texttt{adjG} is specified as the $L \times L$ adjacency matrix, containing only ones and zeros.
    %
    % L is the number of nodes
    % K is the number of external variables
    % order is the polynomial order of the model
    % 
    % topo needs to have 3 fields, G, H, and R.    
    % topo.G is the structure of G, containing a 1 where there is a link and
    % a 0 where there is no link
    % topo.R is the structure of R, containing a 1 where there is a link and
    % a 0 where there is no link
    % topo.H is the structure of H, containing a 1 where there is a link and
    % a 0 where there is no link
    %
    % n needs to have 4 fields, D, NG, NR, and NH
    % n.D is the polynomial order of the D polynomial
    % n.NG is the polynomial order of the N_G polynomial
    % n.NR is the polynomial order of the N_R polynomial
    % n.NH is the polynomial order of the N_H polynomial
    %
    % fixed may contain the three fields G, R, H
    % fixed.G, fixed.R, fixed.H must be specified as either a transfer
    % function matrix, or a real matrix. If specified as a real matrix, then the
    % transfer function with gain as the real matrix will be assigned 

    
    
    
    ni = nargin; % number of input arguments

    % Pre-programmed structures are assigned like this
    if ni == 4 || ni == 5
        L = varargin{2}; % number of nodes
        K = varargin{3}; % number of external inputs
        p = L; % number of noise sources
        
        % Selection of a particular network structure
        switch varargin{1} 
            case 'diagRfullH'
                % Assign orders
                orders.D = varargin{4};
                orders.NG = varargin{4};
                orders.NR = varargin{4};
                orders.NH = varargin{4};
                
                % Assign topology
                topo.NG = ones(L) - eye(L);
                topo.NR = eye(L,K);
                topo.NH = ones(L);
                
                % Make it
                net = network_polynomial(L,K,orders,topo);
                
            case 'diagR'
                % Assign orders
                orders.D = varargin{4};
                orders.NG = varargin{4};
                orders.NR = varargin{4};
                orders.NH = varargin{4};
                
                % Assign topology
                topo.NG = ones(L) - eye(L);
                topo.NR = eye(L,K);
                topo.NH = eye(L);       
                
            case 'fixedR'
                % Assign orders
                orders.D = varargin{4};
                orders.NG = varargin{4};
                orders.NR = varargin{4}
                orders.NH = varargin{4};
                
                % Assign topology
                topo.NG = ones(L) - eye(L);
                topo.NR = eye(L,K);
                topo.NH = eye(L);   
                
                % Fix values of R to 1
                fixed.R = tf(ss([],[],[],eye(L,K),-1));
                
            case 'fixedRfullH'
                % Assign orders
                orders.D = varargin{4};
                orders.NG = varargin{4};
                orders.NR = varargin{4}
                orders.NH = varargin{4};
                
                % Assign topology
                topo.NG = ones(L) - eye(L);
                topo.NR = eye(L,K);
                topo.NH = ones(L);   
                
                % Fix values of R to 1
                fixed.R = tf(ss([],[],[],eye(L,K),-1));
                
                
        end        
        
        % If a custom topology for G was entered, then assign it
        if ni == 5
            topo.NG = varargin{5};
        end
        
    end
    
    % Your own structure 
    if   (2 <= ni) && (ni <= 3)
        
        % Assign orders
        orders = varargin{2};
        
        % Assign topology
        topo = varargin{1};
        
        % Assign variables
        L = size(topo.NG,1); % number of nodes
        K = size(topo.NR,2); % number of external inputs   
        p = size(topo.NH,2); % number of noise sources 
        
        % Assign fixed transfer functions
        if ni == 3
            fixed = varargin{3};
        end
    end

    % Check fixed dynamics 
    % Assign 0 if there are none
    % Convert to transfer function if boolean
    fixed.sure = 1; % This is to create the variable fixed if it doesnt exist
    if ~isfield(fixed,'G')
        fixed.G = tf(ss([],[],[],zeros(L),-1));
    elseif ~isa(fixed.G,'tf')
        fixed.G = tf(ss([],[],[],fixed.G,-1));
    end
    if ~isfield(fixed,'R')
        fixed.R = tf(ss([],[],[],zeros(L,K),-1));
    elseif ~isa(fixed.R,'tf')
        fixed.R = tf(ss([],[],[],fixed.R,-1));
    end
    if ~isfield(fixed,'H')
        fixed.H = tf(ss([],[],[],zeros(L,p),-1));
    elseif ~isa(fixed.H,'tf')
        fixed.H = tf(ss([],[],[],fixed.H,-1));
    end
                    
                    
                    
    
    % This creates a matrix with a 1 at each non-zero transfer function
    fixed_topo.G = boolean(squeeze(sum(impulse(fixed.G).^2,1)));  
    fixed_topo.R = boolean(squeeze(sum(impulse(fixed.R).^2,1)));  
    fixed_topo.H = boolean(squeeze(sum(impulse(fixed.H).^2,1)));  
    
    % Select only modules that contain parameters, i.e. no fixed ones
    parameterized.NG = topo.NG - fixed_topo.G;
    parameterized.NR = topo.NR - fixed_topo.R;
    parameterized.NH = topo.NH - fixed_topo.H;
    
    
    % Make it
    net = network_polynomial(L,K,p,orders,parameterized);

    % Add adjacency matrix
    net.adjacencyG = topo.NG;
    net.adjacencyR = topo.NR;
    net.adjacencyH = topo.NH;

    % Add fixed dynamics
    net.fixedG = fixed.G;
    net.fixedR = fixed.R;
    net.fixedH = fixed.H;
    
    % Add G,R,H,Lambda
    [net.G,net.R,net.H] = model_to_transfer(net,[],net);
    net.Lambda = eye(p);
    net.Gamma = zeros(L-p,p); % Actually the H is H-[0;Gamma]
    
    
    % Test generic network identifiability
    net1.G = net.adjacencyG;
    net1.R = net.adjacencyR;
    net1.H = net.adjacencyH;
    [id,id_n,id_m] = test_identifiability(net1,fixed_topo);
    if id==0
        ind = find(id_n == 0);
        disp(['Warning: Rows ' mat2str(ind) ' of the network model are not generically network identifiable']);
    end
    
end


