function [ nid ] = isotopicCalculator( M, plot )
% 
% calculate isotopic distributions of molecules using the FFT
% http://www.ms-utils.org/isotop.html
% (c) Magnus Palmblad, 1999.  Adapted by James McKenzie 2012, 2014.
%
% M=[378 234 65 75 6]; % empirical formula, e.g. bovine insulin
% fast radix-2 fast-Fourier transform algorithm is used

MAX_ELEMENTS=18; MAX_MASS=2^13;

% the original algorithm has an empty first column.  This i found to be not
% amenable to expansion of the nuclide matrix.  So i took it out, and it
% gets put back in later on!  Right about here...
[A,B] = periodicTable(MAX_ELEMENTS,MAX_MASS);

% Here perform the calculation
[mz,nid] = calcDist(M,A,B,MAX_ELEMENTS,MAX_MASS);

if nargin == 2
    figure; stem(mz,nid);
end    

return

% Enforce output of 5 arguments
if length(nid) > 5
    nid = nid(1:5);
elseif length(nid) < 5
    tmp = zeros(1,5);
    tmp(1:length(nid)) = nid;
    nid = tmp;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B] = periodicTable(MAX_ELEMENTS,MAX_MASS)

% isotopic abundancies stored in A
A=zeros(MAX_ELEMENTS,MAX_MASS-1);
B=zeros(MAX_ELEMENTS,MAX_MASS-1);

% Define the compositions...
A(2,1:2)    =[0.999885 0.000115];                                       % H 
A(1,12:13)  =[0.9893 0.0107];                                           % C
A(4,14:15)  =[0.99636 0.00364];                                         % N
A(5,16:18)  =[0.99757 0.00038 0.00205];                                 % O
A(11,32:36) =[0.9499 0.0075 0.0425 0 0.0001];                           % S
A(3,6:7)    =[0.0759 0.9241];                                           % Li
A(7,23)     =[1];                                                       % Na
A(8,24:26)  =[0.7899 0.100 0.1101];                                     % Mg
A(9,28:30)  =[0.92223 0.04685 0.03092];                                 % Si
A(10,31)    =[1];                                                       % P
A(12,35:37) =[0.7576 0 0.2424];                                         % Cl
A(13,39:41) =[0.932581 0.000117 0.067302];                              % K
A(14,40:48) =[0.96941 0 0.00647 0.00135 0.02086 0 0.00004 0 0.00187];   % Ca
A(15,75)    =[1];                                                       % As
A(16,74:82) =[0.0089 0 0.0937 0.0763 0.2377 0 0.4961 0 0.0873];         % Se
A(17,79:81) =[0.5069 0 0.4931];                                         % Br
A(18,127)   =[1];                                                       % I
A(6,19)     =[1];                                                       % F

% Define the accurate masses...
B(2,1:2)    =[1.0078 2.0141];                                           % H 
B(1,12:13)  =[12.0000 13.0034];                                         % C
B(4,14:15)  =[14.0031 15.0001];                                         % N
B(5,16:18)  =[15.9949 16.9991 17.9991];                                 % O
B(11,32:36) =[31.9721 32.9715 33.9679 0 35.9670];                       % S
B(3,6:7)    =[6.0151 7.0160];                                           % Li
B(7,23)     =[22.9898];                                                 % Na
B(8,24:26)  =[23.9850 24.9855 25.9826];                                 % Mg
B(9,28:30)  =[27.9769 28.9765 29.9738];                                 % Si
B(10,31)    =[30.9738];                                                 % P
B(12,35:37) =[34.9689 0 36.9659];                                       % Cl
B(13,39:41) =[38.9637 39.9640 40.9618];                                 % K
B(14,40:48) =[39.9626 0 41.9586 42.9588 43.9555 0 45.9537 0 47.9525];   % Ca
B(15,75)    =[74.9216];                                                 % As
B(16,74:82) =[73.9225 0 75.9192 76.9199 77.9173 0 79.9165 0 81.9167];   % Se
B(17,79:81) =[78.9183 0 80.9162];                                       % Br
B(18,127)   =[126.9045];                                                % I
B(6,19)     =[18.9984];                                                 % F

% The original algorithm has an empty first column.  This i found to be not
% amenable to expansion of the nuclide matrix.  So i took it out, and it
% gets put back in later on!  Right about here...
A = [zeros(MAX_ELEMENTS,1) A];
B = [zeros(MAX_ELEMENTS,1) B];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mass,nid] = calcDist(M,A,B,MAX_ELEMENTS,MAX_MASS)

% Trying to calculate the mass distributions too, but to no avail.
[idA] = coreCalc(M,A,MAX_ELEMENTS,MAX_MASS);
%[idB] = coreCalc(M,B,MAX_ELEMENTS,MAX_MASS)

% Transform to 100 for max intensity
mx = max(idA);
nid = 100 .* idA ./ mx;

% Remove peaks at less than 0.1% base peak intensity
[x] = find(nid > 0.1);
nid = nid(x);

% Put the distribution into a mono-spaced vector
mass = x - x(1);
numP = max(mass) + 1;
dist = zeros(1,numP);
for n = 1:numel(nid)
    
    ind = mass(n) + 1;
    dist(1,ind) = nid(n);
    
end
   
nid = dist;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [id] = coreCalc(M,A,MAX_ELEMENTS,MAX_MASS)
% These are done for both intensity and mass

% FFT along each element's isotopic distribution
tA=fft(A,[],2);
%tA = A;

% multiply transforms (elementwise)
ptA = ones(1,MAX_MASS);
for i = 1:MAX_ELEMENTS,
  ptA = ptA .* (tA(i,:) .^ M(i));         
end

% Inverse FFT to get convolutions
riptA = real(ifft(ptA));

% Shift to real mass
id=zeros(1,MAX_MASS);
id(1:MAX_MASS-1) = riptA(2:MAX_MASS);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%