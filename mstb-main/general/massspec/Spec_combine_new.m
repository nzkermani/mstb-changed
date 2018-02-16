function [MZ,Output]=Spec_combine(OrbData,resolution,min_mz,max_mz)
% Takes output from imzML reading script 'Data_extract' and combines all
% spectra from one sample

%Input is OrbData cell array
no_samples=size(OrbData,2);
MZ=min_mz:resolution:max_mz;

Output=zeros(no_samples,length(MZ));

% Define noise regions...
FTnoiseregions=[197.127000000000;350.717000000000;398.801000000000;...
    477.365000000000;533.774000000000;533.802000000000;...
    542.812000000000;687.417000000000;781.647000000000;...
    792.115000000000;1842.69100000000;1909.46700000000;...
    1984.74200000000;1999.56000000000];
%FTn=round(FTnoiseregions);
%FTnn(:,1)=FTn-5;
%FTnn(:,2)=FTn+5;

%no_samples=2;
%n=2;
%F=1;
%TIC2=1;


% Holding={};

% each of the analyses is one sample of x spectra process each one by one
% then put them on a combined m/z scale in the final stage
for i=1:no_samples;
    
    % Get the orb data for a single scan...
    Sample = OrbData{i};
    
    % How many spectra?
    no_spectra = size(Sample,2);
    
    
    Good=1;
    TIC2=1;
    
    
    for l=1:no_spectra;
        X=Sample{l};
        T=X(:,2);Tic=sum(T);
        TIC2(l,:)=Tic;
    end
    
    med1=median(TIC2);
    Numbers=1:1:no_spectra;
    G2=TIC2>(0.5*med1);
    Good=Numbers(G2);
    
    Temp=zeros(length(Good),length(MZ));
    
    
    
    
    for j=1:length(Good);
        
        l=Good(j);
        X=Sample{l};
        
        mask = X(:,1) >= min_mz & X(:,1) <= max_mz;
        X = X(mask,:);       
                
        
        % Original interpolation line
        G1 = interp1(X(:,1),X(:,2),MZ);
        
        %G2 = interp1(X(:,1),X(:,2),MZ,'nearest');
        
        
%         figure; hold on; 
%         stem(X(:,1),X(:,2),'b'); 
%         stem(MZ,-G1,'r');
%         stem(MZ, G2,'m');
        
%         % Binning step        
%         MZ0=X(:,1)';
%         Spectra=X(:,2)';        
%         nospec=size(Spectra,1);
%         G1=zeros(nospec,length(MZ));        
%         Array=zeros(1,length(MZ0));
%         for k=1:length(MZ0);            
%             a=MZ0(k);            
%             H=(MZ-a).^2;            
%             [C1,C2]=min(H);            
%             Array(:,k)=C2;
%         end
%         
%         
%         
%         for m=1:length(MZ0);
%             
%             
%             G1(1:nospec,Array(1,m))=G1(1:nospec,Array(1,m))+Spectra(1:nospec,m);
%             
%         end
        
        
        %G1=interp1(X(:,1),X(:,2),MZ);
        
        Temp(j,:)=G1';
        
    end
    
    Output(i,:)=mean(Temp);
    
    
    
end

Output(isnan(Output))=0;
Output(isinf(Output))=0;

end









