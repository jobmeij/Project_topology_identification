function [identifiable,identifiable_nodes,identifiable_modules] = test_identifiability(varargin)
% Function to test the generic network identifiability of the network
% model. This is tested by checking conditions on the maximum number of 
% vertex disjoint paths. 
%
% A way to check the maximum number of vertex disjoint paths is splitting 
% nodes into an input and an output node, such that it equals the maximum 
% number of edge disjoint paths. The maximum number of edge disjoint paths 
% is the maximum flow of the unweighted graph.
%
% The function can be used in this way
% [ide,id_nodes,id_modules] = test_identifiability(adjacency,fixed)
% adjacency has the fields .G, .R, .H, which each contain the adjacency matrix
% of G, R and H
% fixed is optional, and has the optional fields .G, .R, .H, which each
% specify which modules of G, R and H are non-zero but fixed/known

    % Shortcuts
    adjacency = varargin{1}; % adjacency matrices
    L = size(adjacency.G,1); % number of nodes
    K = size(adjacency.R,2); % number of external excitations
    p = size(adjacency.H,2); % rank of noise model
    
    % these are non-zero transfers that do not contain parameters in the model
    fixed_dynamics.G = 0 * adjacency.G;
    fixed_dynamics.R = 0 * adjacency.R;
    fixed_dynamics.H = 0 * adjacency.H;   
    
        
    if nargin > 1 && isfield(varargin{2},'G')
        fixed_dynamics.G = varargin{2}.G; 
    end 
    if nargin > 1 && isfield(varargin{2},'R')
        fixed_dynamics.R = varargin{2}.R; 
    end 
    if nargin > 1 && isfield(varargin{2},'H')
        fixed_dynamics.H = varargin{2}.H; 
    end
    
    
    
    
    %% Build equivalent network
    % Build a network where r and e are nodes, and the nodes w are split
    % into "w_in" and "w_out": the nodes are [w_in; w_out; r; e; s; t]
    
    % Some shortcuts to label matrix indices
    Ltotal = 3*L + K + 2;
    ind_in = 1:L;
    ind_out = L+1:2*L;
    ind_r = 2*L+1:2*L+K;
    ind_e = 2*L+K+1:2*L+K+p;
    ind_s = 2*L+K+p+1;
    ind_t = 2*L+K+p+2;
    
    Atotal = zeros(Ltotal,Ltotal); % empty adjacency matrix
    
    Atotal(ind_out,ind_in) = eye(L); % Connect w_in to w_out because the nodes were split
    
    Atotal(ind_in,ind_out) = adjacency.G; % Connect w_out to w_in by the original adjacency matrix
    
    Atotal(ind_in,ind_r) = adjacency.R; % Connect r to w_out by the original adjacency matrix
    Atotal(ind_in,ind_e) = adjacency.H; % Connect e w_out by the original adjacency matrix
    
    
    %G = digraph(Atotal',{'w1in','w2in','w3in','w4in','w5in','w1out','w2out','w3out','w4out','w5out','r1','r2','r3','r4','r5','e1','e2','e3','e4','e5','s','t'});
    %G = digraph(Atotal',{'w1in','w2in','w3in','w1out','w2out','w3out','r1','r2','e1','e2','e3','s','t'});
    %plot(G)
    
    
    
    %% Test identifiability of each row
    for node = 1:1:L
        
        % Build adjacency matrix to test identifiability for row j 
        Ast = Atotal;
        
        
        % Select set of excitations for which R_{jk} does not contain
        % parameters, and for which H_{jl} does not contain parameters
        fix = [];
        if size(fixed_dynamics.R,1) > 0
            fix = [fix fixed_dynamics.R(node,:)];
        end
        if size(fixed_dynamics.H,1) > 0
            fix = [fix fixed_dynamics.H(node,:)];
        end
        U = Atotal(node,[ind_r ind_e]) == 0 + logical(fix);
        Ast([ind_r ind_e],ind_s) = U';
        
        
        % Select set of k for which G_{jk} contains parameters
        Y = Atotal(node,ind_out) == 1 - logical(fixed_dynamics.G(node,:));
        Ast(ind_t,ind_out) = Y;
        
        %G = digraph(Ast',{'w1in','w2in','w3in','w4in','w5in','w1out','w2out','w3out','w4out','w5out','r1','r2','r3','r4','r5','e1','e2','e3','e4','e5','s','t'});
        %plot(G)
        
        % max flow = edge disjoint paths from s to t = vertex disjoint
        % paths from U to Y
        Gst = digraph(Ast');
        mf = maxflow(Gst,ind_s,ind_t);
        
        % If sufficient paths, then node is identifiable
        identifiable_nodes(node) = mf >= sum(Y);
        
        if identifiable_nodes(node)
            identifiable_modules(node,:) = ones(1,L);
        else
            for ii = 1:1:L
                if Y(ii)
                    Ybar = Y;
                    Ybar(ii) = 0;
                    
                    Ast_bar = Ast;
                    Ast_bar(ind_t,ind_out) = Ybar;
                    Gst_bar = digraph(Ast_bar');
                    mf_bar = maxflow(Gst_bar,ind_s,ind_t);
                    
                    identifiable_modules(node,ii) = mf > mf_bar;
                else
                    identifiable_modules(node,ii) = 1;
                end
            end
        end
        
    end
    
    % If all nodes identifiable, then network is identifiable
    identifiable = sum(identifiable_nodes) == L;
end