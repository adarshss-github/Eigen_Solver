function [Phi,f] = subit(K,M,p,tol)

%Function for subspace iteration (Will extract the first p undamped modes of a large dof dynamic system)
%By Adarsh S (PhD Candidate IIT Kanpur)

%Example [Phi,f] = subit(K,M,7,1e-3)

%Inputs and Outputs of the function
% K - Stiffness Matrix
% M - Mass Matrix
%p - No. of Modes to be extracted
%tol - Tolerence for convergence                                                            
%Phi - First p mode shapes (Mass ortho-normalized)
%f - First p natural frequencies in Hz
% The eigen values and vectors for the smaller eigen value problem is found by using the inbuilt matlab function which seems to be stpuid,
%but this program was originally inteded as a pedagogical tool

 %%
 
 LK = chol(K,'lower') ; %Cholesky decomposition of stiffness matrix
 [n1,~] = size(K) ;
 [n2,~] = size(M) ;
                              
 if n1~=n2
     fprintf(' Error: Dimensions of K and M matrices not matching \n' ) ;
 end
 
 m1 = 2*p ;     %Criteria for number of starting vectors as given by Bathe
 m2 = p + 8 ;
 
 if m1>m2
     m = m2 ;
 else
     m = m1 ;
 end
 
 if m>n1
     m = n1 ;
 end
     
  ksmall = diag(K) ;
  msmall = diag(M) ;
  
  r = ksmall./msmall ;
  [~,I] = sort(r) ;
  
  Vi = zeros(n1,m) ;
  Vi(:,1) = msmall ;
  dum = zeros(n1,1) ;
  
  for i = 2:1:m                  %Formation of initial matrix for iteration as given by Bathe
      dum(I(i-1),1) = 1 ;
      Vi(:,i) = dum ;
      dum = zeros(n1,1) ;
  end

  %%
  %Start of initial step of iteration
  V1 = zeros(n1,m) ;
  V2 = zeros(n1,m) ;
  
  for i = 1:1:m                 %Inverse iteration step
       V1(:,i) = LK'\(LK\Vi(:,i)) ;
  end
  
  d1 = zeros(n1,1) ;
  A = zeros(m,m) ;
  B = zeros(m,m) ;
  A1 = zeros(m,m) ;
  Vmn = zeros(m,m) ;
  Kstar = zeros(m,m) ;
  Mstar = zeros(m,m) ;
  c = 0 ;
  k = 0 ;
  
  Kstar = V1'*K*V1 ;
  Mstar = V1'*M*V1 ;
  
  [A,B] = eig(Kstar,Mstar) ;   % Subspace eigen value problem
  d1 = diag(B) ;
  [d1,I] = sort(d1) ;
  
  for i = 1:1:m
      A1(:,i) = A(:,I(i)) ;
  end
  
   for i = 1:1:m                 %Mass ortho-normalization
      c  = sqrt(A1(:,i)'*Mstar*A1(:,i)) ;
      Vmn(:,i) =  A1(:,i)/c ;
   end
 
  V2 = V1*Vmn ;
  j = 0 ;
  
 %Start of main iteration
 while j<p
     
 j = 0 ;
 
 for i = 1:1:m             %Inverse iteration step
   V1(:,i) = LK'\(LK\V2(:,i)) ;
 end 
  
 Kstar = V1'*K*V1 ;
 Mstar = V1'*M*V1 ;
  
 [A,B] = eig(Kstar,Mstar) ;  % Subspace eigen value problem
 d2 = diag(B) ;
 [d2,I] = sort(d2) ;
 
  for i = 1:1:m
      A1(:,i) = A(:,I(i)) ;
  end
  
   for i = 1:1:m                     %Mass ortho-normalization
      c  = sqrt(A1(:,i)'*Mstar*A1(:,i)) ;
      Vmn(:,i) = A1(:,i)/c ;
   end
  
   V2 = V1*Vmn ;
   
   for i = 1:1:p                     %Convergence check based on tolerence given
       k = abs((d2(i)-d1(i))/d2(i)) ;
       if k<tol
           j = j+1 ;
       else
           break
       end
   end
   
   d1 = d2 ;
   
 end
 
Phi = V2(:,1:p) ;
d = d2(1:p) ;
f = sqrt(d)/(2*pi) ;

end
  
  
  