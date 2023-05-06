function H=HessMp(f,x0,varargin)
% computes the Hessian matrix of f evaluated at x0. If x0 has K elements,
% function returns a KxK matrix. The function f, given the ML context, is
% expected to be a column vector
% uses central differences 

% f should return either a scalar or a column vector
% x0 should be a column vector of parameters
f0=feval(f,x0,varargin{:}); 
[T,co]=size(f0);
if co>1; error('Error in HessMp, The function should be a column vector or a scalar'); end

[k,c]=size(x0);
if k<c,
    x0=x0';
end
k=size(x0,1); % number of parameters wrt which one should compute gradient

h=0.00001; %some small number

H=zeros(k,k); %will contain the Hessian
e=eye(k); 

h2=h/2;
for ii=1:k;
      if x0(ii)>100; % if argument is big enough, compute relative number   
        x0P= x0.*( ones(k,1) +  e(:,ii) *h2 );
        x0N= x0.*( ones(k,1) -  e(:,ii) *h2 );
        Deltaii = x0(ii)*h;
    else
        x0P = x0 +  e(:,ii) *h2;
        x0N = x0 -  e(:,ii) *h2;
        Deltaii = h;
    end
    
    for jj=1:ii
    if x0(jj)>100; % if argument is big enough, compute relative number   
        x0PP = x0P .* ( ones(k,1) +  e(:,jj) *h2 );
        x0PN = x0P .* ( ones(k,1) -  e(:,jj) *h2 );
        x0NP = x0N .* ( ones(k,1) +  e(:,jj) *h2 );
        x0NN = x0N .* ( ones(k,1) -  e(:,jj) *h2 );
        Delta = Deltaii*x0(jj)*h;
    else
        x0PP = x0P  +  e(:,jj) *h2; 
        x0PN = x0P  -  e(:,jj) *h2; 
        x0NP = x0N  +  e(:,jj) *h2; 
        x0NN = x0N  -  e(:,jj) *h2; 
        Delta = Deltaii*h;
    end
    
        fPP = feval(f,x0PP,varargin{:});   % forward,forward
        fPN = feval(f,x0PN,varargin{:});   % forward,backward
        fNP = feval(f,x0NP,varargin{:});    % backward,forward
        fNN = feval(f,x0NN,varargin{:});    % backward,backward
        
        H(ii,jj)=(sum(fPP)-sum(fPN)-sum(fNP)+sum(fNN))/Delta;
        H(jj,ii)=H(ii,jj);
    end
end