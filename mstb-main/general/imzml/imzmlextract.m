function [MZ,Int,info] = imzmlextract(binning_resolution,lowmass,highmass)
% Gets the required information from an imzML file such that the processed data 
% in the ibd file can be harvested using imzML_rip. 

% Emrys Jones 

resolution=binning_resolution;
disp('Select .imzML file');


[filename,filepath] = uigetfile('*.imzML');

if ~strcmp(filepath(end),filesep)
    filepath = [filepath filesep];
end

filename = [filepath filename];

um=lowmass:resolution:highmass;


fid=fopen(filename);

C = textscan(fid, '%s',  'delimiter', '##$' , 'multipledelimsasone', 1);



D=C{:};



clear C 
f1='<referenceableParamGroup id="mzArray">';
f2='<referenceableParamGroup id="intensityArray">';

for i=1:length(D);
    if strcmp(f1,D(i))==1;
        in1=i;

    elseif strcmp(f2,D(i))==1;
        in2=i;
    end
end

%in1=in1+2;
%in2=in2+2;
in1=in1+4;
in2=in2+4;

mzAprec=D{in1};mzAprec(length(mzAprec)-11:length(mzAprec))=[];mzAprec(1:49)=[];
intAprec=D{in2};intAprec(length(intAprec)-11:length(intAprec))=[];intAprec(1:49)=[];


if mzAprec=='64-bit float';
    preMZ='float64';
    
else if mzAprec=='32-bit float';
        preMZ='float32';
        
    else preMZ='double';
        
    end
end

if intAprec=='64-bit float';
    preINT='float64';
    
else if intAprec=='32-bit float';
        preINT='float32';
        
    else preINT='double';
        
    end
end

if intAprec=='64-bit float';
    preINT='float64';
    
else if intAprec=='32-bit float';
        preINT='float32';
        
    else preINT='double';
        
    end
end



% number of pixels

%XP='<cvParam cvRef="IMS" accession="IMS:1000042" name="max count of pixel x" value="104"/>';
%YP='<cvParam cvRef="IMS" accession="IMS:1000043" name="max count of pixel y" value="139"/>';

XP='<cvParam cvRef="IMS" accession="IMS:1000042" name="max count of pixel x" value="104"/>';
YP='<cvParam cvRef="IMS" accession="IMS:1000043" name="max count of pixel y" value="139"/>';

 
for i=1:length(D);
    if strncmpi(XP,D(i),79)==1;
        pixelsX=i;
        break
    end
end


    
    



pixelsX=D{pixelsX};

maxd=size(pixelsX,2);
pixelsX(maxd-2:maxd)=[];
pixelsX(1:81)=[];
NoXP=str2num(pixelsX);



for i=1:length(D);
if strncmpi(YP,D(i),79)==1;
pixelsY=i;
break 
end
end

pixelsY=D{pixelsY};

maxd=size(pixelsY,2);
pixelsY(maxd-2:maxd)=[];
pixelsY(1:81)=[];
NoYP=str2num(pixelsY);


%find the first pixel

px1='<cvParam cvRef="IMS" accession="IMS:1000050" name="position x" value="1"/>';

py1='<cvParam cvRef="IMS" accession="IMS:1000051" name="position y" value="1"/>';

ps1='<cvParam cvRef="IMS" accession="IMS:1000103" name="external array length" value="289505"/>';

for i=1:length(D);
if strcmp(px1,D(i))==1;
PX1=i;
break 
end
end




for i=1:length(D);
if strcmp(py1,D(i))==1;
PY1=i;
break 
end
end

for i=1:length(D);
if strncmpi(ps1,D(i),78)==1;
PS1=i;
break 
end
end


px2='<cvParam cvRef="IMS" accession="IMS:1000050" name="position x" value="2"/>';

for i=1:length(D);
if strcmp(px2,D(i))==1;
PX2=i;
break 
end
end

Shift1=PX2-PX1;

for i=1:length(D);
if strncmpi(px1,D(i),70)==1;
End=i;

end
end

NoS=End-PX1;

Spectra=1+(NoS/Shift1);


%Spectra=NoXP*(NoYP-1); % added in -1 for specific data fix

PixX=zeros(Spectra,1);
PixY=zeros(Spectra,1);
SpecL=zeros(Spectra,1);



for i=1:Spectra;
    
P_X(i,:)=D((PX1+(Shift1*(i-1))),:);         % normally plus 27
P_Y(i,:)=D((PY1+(Shift1*(i-1))),:); 
SL(i,:)=D((PS1+(Shift1*(i-1))),:); %(normally plus 6)



end

for i=1:Spectra;

YP=cell2mat(P_Y(i,:)); maxc=size(YP,2);YP(maxc-2:maxc)=[];YP(1:70)=[];PixY(i,:)=str2num(YP);

XP=cell2mat(P_X(i,:)); maxc=size(XP,2);XP(maxc-2:maxc)=[];XP(1:70)=[];PixX(i,:)=str2num(XP);

SLP=cell2mat(SL(i,:)); maxc=size(SLP,2);SLP(maxc-2:maxc)=[];SLP(1:81)=[];SpecL(i,:)=str2num(SLP);

end


% This extracted data used to populate a datacube from the ibdfile. 


filename(length(filename)-4:length(filename))=[];
filename=[filename,'ibd'];


fid=fopen(filename,'r'); 


spec=size(SpecL,1);

% masses=zeros(spec,M);
% spectra=zeros(spec,M);


kill=fread(fid,2,'double');
clear kill






Int=zeros(spec,length(um));



for i=1:spec;
            
    D=SpecL(i,:);    
          
    m1=fread(fid,D,preMZ);
    
    X=m1>lowmass & m1<highmass;
    
    m2=m1(X);
    
    
      
    s1=fread(fid,D,preINT);
    
    s2=s1(X);
    
    for j=1:length(s2);
        
        
    [I1,I2]=min((um-m2(j)).^2);
        
       Int(i,I2)=Int(i,I2)+s2(j);
    
       disp(i);
        
    end
   
   %Int(i,:)=interp1(m2,s2,um);
    
    
end


    

    
    
      
  


MZ=um;

info=[NoYP,NoXP];

end



