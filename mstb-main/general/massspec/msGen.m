function [mzVec,spec,opts] = msGen
%msGen - generation of spectra

% Define m/z range and resolution
opts.mzVal = [45 270];
opts.mzRes = 0.01;
opts.fwhm  = 50000;

% Noise options
opts.varPeaks = 0.01;    % ±% variation in peak intensity
opts.noise = 100;     % maximum intensity of noise
opts.posVar = 5;       % ± maximum ppm positional variation of peaks

% How many molecules and spectra to be calculated?
numSpec = 50;

% Read in the molecules: molecular formula, isotopic dists etc...
[mol,ints] = molecules(numSpec);
numMol = numel(mol.form);

% Determine the m/z vector
mzVec = opts.mzVal(1):opts.mzRes:opts.mzVal(2);
numMZ = numel(mzVec);

% Pre-allocate matrix with noise
spec = randn(numSpec,numMZ) * opts.noise;
%spec = spec - min(spec(:));
%spec = zeros(numSpec,numMZ);

% Run through each sample...
for s = 1:numSpec
        
    % Run through each molecule
    for m = 1:numMol
        
        % Multiply each peak by its concentration
        tmp = ints(s,m) * mol.dist{m,1};
        numP = numel(tmp);
        
        % For each peak generate the gaussian distribution, and then add it
        % into the spectrum
        for p = 1:numP
            
            % Modify according to the random variation in peak intensity
            rnd = 1 + (rand(1,numel(tmp)) * (2*0.05)) - 0.05;
            tmp = tmp ./ rnd;
            
            % Generate gaussian
            for g = 1:numel(tmp)
                if tmp(g) > 0
                    
                    % What is the position of this peak?
                    peakMZ = mol.mono(m) + (g-1);
                    
                    % Add in peak uncertainty...
                    posUnc = (2 * rand(1,1)) - 1;                    
                    posUnc = posUnc * opts.posVar * peakMZ / 1e6;                    
                    peakMZ = peakMZ + posUnc;
                                        
                    % Round to fit into current resolution
                    peakMZ = round(peakMZ / opts.mzRes) * opts.mzRes;
                    
                    % The sigma for this peak is fixed by the FWHM and
                    % depends on the peak's m/z
                    sigma = peakMZ / opts.fwhm;
                    
                    % Peak st/fn positions?
                    st = peakMZ - (opts.mzRes * 75);
                    fn = peakMZ + (opts.mzRes * 75);
                    
                    % Calc Gaussian function
                    [y] = genGauss(st,fn,...    % peak start/end mz values
                        opts.mzRes,...          % mz resolution
                        peakMZ,...              % peak position
                        sigma,...               % peak width
                        tmp(g));                % peak intensity
                    
                    % Now determine location in the spectral matrix...
                    stI = 1+round((st-opts.mzVal(1))/opts.mzRes);
                    fnI = 1+round((fn-opts.mzVal(1))/opts.mzRes);
                    
                    % Finally add into the matrix
                    spec(s,stI:fnI) = spec(s,stI:fnI) + y;                    
                    
                end
            end                
        end        
    end    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,ints] = molecules(numSpec)
% List the names / formulae of the molecules

m.name = {'Creatinine','Hippuric_acid',...
    'Citric_acid','Glycine',...
    'TMAO','Histidine',...
    'phenylacetylglutamine','Taurine',...
    'glycolic_acid','Formate'};

m.form = {'C4H7N3O','C9H9NO3',...
    'C6H8O7','C2H5NO2',...
    'C3H9NO','C6H9N3O2'...
    'C13H16N2O4','C2H7NO3S',...
    'C2H4O3','CH2O2'};

numI = numel(m.name);
for n = 1:numI
    
    % Parse the formula
    pf = parseFormula(m.form{n});
    
    % Add a hydrogen ion (ignore mass of electron)
    pf.qty(2) = pf.qty(2) + 1;
    
    % Calculate the monoisotopic mass
    m.mono(n,1) = sum(pf.mono .* pf.qty);
    
    % Calculate isotopic distribution...
    m.dist{n,1} = isotopicCalculatorLite(pf.qty);
    
end

% Intensity matrix...
ints = zeros(numSpec,numI);

% Now read in their concentrations from the text file
fid = fopen('Concentrations_Mix1.txt');
for n = 1:numI
    
    % Read in the molcule name...
    name = textscan(fid,'%s',1);
    
    % Let's match the names...
    cmp = strcmpi(m.name,name{1});
    %fx = find(cmp == 1);
    
    % Now read in the intensities
    tmp = textscan(fid,'%f',numSpec);       
    ints(:,cmp) = tmp{1};
    
end
fclose(fid);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = genGauss(st,fn,tol,pos,sigma,iii)
% Generate a gaussian function...

% Generate the x data...
x = st:tol:fn;

numer = (x - pos) .^ 2;
denom = 2 * (sigma ^ 2);

y = iii * exp(-numer/denom);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%