
% Date : 13/07/2023
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

%% Sub-function of : Main_code.m
%  Sub functions   : shape_fn.m
%                  : lgwt.m

% Input :
%       Finite layer parameters (cp, cs, roS, h)
%       nDivs - Nomber of sub layer in each layer
%       d -  Order of Polynomial
% Output :
%       Stiffness matrices for finite layer (Kxx, Kxz, Kzz, M)


function [A, B, C, M] = stiffness(cs, cp, roS, nElem, h, d)


n     = d + 1;  % Number of Gauss points

% Calculation of n-numbers of integration points and corresponding weights
[GP]  = Lagr_wt(n);  % GP is 2*n matrix; 1st row containg integral point co-ordinate and 2nd row contain weights

% Note : Order of shape function elements (d) - linear(1); quadratic(2);
%        cubic(3); quartic(4);
%        Ordert of stiffness matrix element (order of N'N) (p) - linear(2);
%        quadratic(4); cubic(6); quartic(8);
%        Required integration points p = 2n-1
%        n-point Gaussian quadrature yield an exact result for polynomials of degree 2n − 1 or less

pt    = d + 1;

nL    = length(cs);             % Number of finite layers
NumL  = sum(nElem);             % total number of thin layers
DOF   = 2*pt*NumL-2*(NumL-1);   % Each node having two degree of freedom
A   = zeros(DOF,DOF);         % Allocating space for the stiffness matrices
B   = A; C = A; M = A;

mu    = roS.*cs.^2;             % Lamé's first parameter
Lam   =  roS.*cp.^2-2*mu;       % Lamé's second parameter

% Calculating some index values to perform globalization
e1    = 2*pt -1;
e2    = e1-1;
Ldis  = [1  e1:e2:length(A)];
Ldis1 = Ldis(2:end)+1;

jj = 0;

%%
for j = 1 : nL  % starting with each layer
    
    % Material sub-matrices. Please refer to equation(5) of the paper
    Kxx = [Lam(j)+2*mu(j) 0 ; 0 mu(j)];
    Kxz = [0 Lam(j) ;  mu(j) 0];
    Kzz = [mu(j) 0 ; 0 Lam(j)+2*mu(j)];
    
    for k = 1 : nElem(j) % for each sub layer
        
        jj = jj+1;
        
        for m = 1 : size(GP,2) % at each intrigation point of the sub layer
            
            p_wt = GP(:,m);    % Taking co-ordinate and corrosponding weight of m th integration point
            p  = p_wt(1);      % co-ordinate of the m th integration point
            wt = p_wt(2);      % weight of the m th integration point
            
            [N1,B1] = shape_fn(d,p); % calculating the values of shape funtion (N1) and its derivative (B1)
            % at the integration point p. The size of N1 will be (d+1)*1
            
            N  = double(N1');
            N  = kron(N,eye(2));     % performing Kronecker product as mentioned above equation (12)
            Bs = double(B1');
            L  = h(j)/nElem(j);      % thickness of the layer i
            J  = L/2;
            N_  = Bs/J; N_ = kron(N_,eye(2));
            rc = Ldis(jj):Ldis1(jj);   % index values to perform globalization
            
            %% Calculate elements of equation (11) in paper
            
            A(rc,rc) = A(rc,rc)   + (N.'*Kxx*N) * (J*wt);
            B(rc,rc) = B(rc,rc)   - (N.'*Kxz*N_) * (J*wt);
            B(rc,rc) = B(rc,rc)   + (N_.'*Kxz.'*N)*J*wt;
            C(rc,rc) = C(rc,rc)   + (N_.'*Kzz*N_) * (J*wt);
            M(rc,rc)   = M(rc,rc) + (N.'*eye(2)*N) * (J*wt*roS(j));
        end
    end
end


end

%% Sub-function of : Stiffness.m
% Computes the Legendre-Gauss nodes (m) and weights (w)  on an 
% interval [+1 -1]

% Input :
%       n - Gauss points

% Output :
%       GP - integration points and corresponding weights


function [GP] = Lagr_wt(n)

n  = n-1;
n1 = n+1; n2 = n+2;
V  = zeros(n1,n2); % Vandermonde Matrix
Vd = zeros(n1,n2); % Derivative
xu = linspace(-1,1,n1)';
m  = cos((2*(0:n)'+1)*pi/(2*n+2)) + (0.27/n1) * sin(pi*xu*n/n2); % initial uess
m1 = 2.*ones(size(m));

while max(abs(m-m1)) > eps    
    V(:,1) = 1;    V(:,2) = m;  
    for k = 2:n1
        V(:,k+1) = ((2*k-1)*m.*V(:,k)-(k-1)*V(:,k-1)) / k;
    end    
    Vd = n2 * (V(:,n1)-m.*V(:,n2))./(1-m.^2);
    m1 = m;
    m  = m1-V(:,n2)./Vd;    
end
w  = 2./((1-m.^2).*Vd.^2)*(n2/n1)^2;
m  = flip(m'); w = flip(w');
GP = [m; w];

end
%%