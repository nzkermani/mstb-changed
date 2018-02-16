%% Alignment algorithm
% Find the smallest m/z in the samples place them in the array A
topTheQueue = mz1;
% Find the smallest in list A, Ax 
% Form a new bin BA place Ax in it
A(1) = min(min(topTheQueue));
% Replace Ax from the according sample in array A
[i, j] = find(topTheQueue == A(1));
topTheQueue(i,j) = mz(i,j,2);
% Add m/zs that fall 2*Eppm from Ax
Eppm = 0.000008*A(1);
temp = topTheQueue(topTheQueue-A(1) < 2*Eppm & topTheQueue-A(1) == 2*Eppm );
% Replace m/zs from the according sample in array A
if(isempty(temp))
    test = 1;
end

% Find the smallest m/z in the samples place them in the array B
% Find the smallest in list B, Bx 
% Form a new bin BB place Bx in it
B(1) = min(min(topTheQueue));
% Replace Bx from the according sample in array B
[i, j] = find(topTheQueue == B(1));
topTheQueue(i,j) = mz(i,j,2);
% Add m/zs that fall 2*Eppm from Bx
Eppm = 0.000008*B(1);
temp = topTheQueue(topTheQueue-B(1) < 2*Eppm & topTheQueue-B(1) == 2*Eppm );
% Replace m/zs from the according sample in array B
if(isempty(temp))
    test = 1;
end

% Check if the largest m/z value value of BA and smallest m/z alue of BB
% are within " * Eppm of each other
if(abs(A(1)-B(1))<2*Eppm  | abs(A(1)-B(1))==2*Eppm)
    test = 1;
else
   % bins are noit within range 
   % check if bins BA has peaks that are more than 2*Eppm from the mean of
   % the bin and reassign m/z peaks to BA and BB using 4i-4v
end

output.mz = A(1);
A(1) = B(1);

% repeat steps 3,4,5,6 until no more data left in the samples
    
