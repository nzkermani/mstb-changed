function [ output_args ] = imzmlJSON(fp,fn,varargin)
%imzmlJSON - create a JSON file for any imzML file so that it can be
%uploaded to the METASPACE engine using the batch functionality.

% Get the stuff from the input arguments
[opts] = readArgsData(varargin);

% Now we can print out the JSON file...
[json] = jsonParts(fn,opts);

fid = fopen([fp fn '.json'],'w');
fprintf(fid,'%s',json);
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = readArgsData(argsin)
% Read the arguments and then the data if it wasn't passed

% Define the defaults here
opts.organism   = 'Homo sapiens (human)';
opts.orgPart    = '';
opts.sampStab   = 'Fresh frozen';
opts.condition  = 'Diseased';
opts.ionSource  = 'DESI';
opts.analyser   = 'Orbitrap';
opts.database   = '"HMDB", "ChEBI"';
opts.polarity   = '';
opts.adducts    = [];
opts.rpMZ       = 200;
opts.rp         = 100000;
opts.sampDesc   = '';
opts.addInfo    = '';


% Run through each pair
nArgs = length(argsin);
for i = 1:2:nArgs
    
    switch lower(argsin{i})
        
        case 'organism'
            tmp = lower(argsin{i+1});
            switch tmp
                case 'homo sapiens (human)'
                otherwise
                    error('Only human samples expected');
            end
            
        case {'part','organismpart','organism part'}
            tmp = argsin{i+1};
            switch lower(tmp)
                case 'brain'
                    opts.orgPart = 'Brain';
                case 'liver'
                    opts.orgPart = 'Liver';
                case {'colon','colorectal'}
                    opts.orgPart = 'Colon';
                case 'breast'
                    opts.orgPart = 'Breast';
                case {'ovarian','ovary'}
                    opts.orgPart = 'Ovary';
                case {'lymph','lymphnode','lymph node'}
                    opts.orgPart = 'Lymphnode';
                otherwise
                    error('Unknown organism part');
            end
            
        case 'condition'
            opts.condition = argsin{i+1};
            
        case {'stabilisation','samplestab','sample stabilisation'}
            tmp = lower(argsin{i+1});
            switch tmp
                case 'fresh frozen'
                    opts.sampStab = 'Fresh frozen';
                otherwise
                    error('Unknown stabilisation type')
            end
                        
        case {'source','ionisation','ionsource','ion source'}
            tmp = argsin{i+1};
            switch lower(tmp)
                case 'desi'
                    opts.ionSource = 'DESI';
                case 'maldi'
                    opts.ionSource = 'MALDI';
                otherwise
                    error('Unfamiliar ionisation source');
            end
            
        case 'analyser'
            tmp = argsin{i+1};
            switch lower(tmp)
                case 'orbitrap'
                    opts.analyser = 'Orbitrap';                    
                case 'tof reflector'
                    opts.analyer = 'TOF reflector';                    
                otherwise
                    error('Unfamiliar analyser');
            end
            
        case {'rpmz','mz','rp mz'}
            if isnumeric(argsin{i+1})
                opts.rpMZ = argsin{i+1};
            else
                error('Bad mz for resolving power specified');
            end
            
        case {'rp','resolving power','resolvingpower'}
            if isnumeric(argsin{i+1})
                opts.rp = argsin{i+1};
            else
                error('Bad resolving power specified');
            end
            
        case 'polarity'
            tmp = argsin{i+1};
            switch lower(tmp(1))
                case 'p'
                    opts.polarity = 'Positive';
                case 'n'
                    opts.polarity = 'Negative';
                otherwise
                    error('Unknown polarity');
            end
            
        case 'adducts'
            % Just assume that these are correct...
            opts.adducts = argsin{i+1};
            
        case 'database'
            disp('Default to using HMDB');
            
        case {'sample description','description'}
            opts.sampDesc = argsin{i+1};
            
        case {'additional information','information','extra'}
            opts.addInfo = argsin{i+1};
            
        otherwise
            disp(argsin{i})
            error('Entry not recognised - please review');
    end    
end

% Must specify an organism part
if isempty(opts.orgPart)
    error('Organism part not specified');
end

% Check that polarity is specified correctly
if isempty(opts.polarity)
    error('Polarity not specified');
end

% Need to add the default adducts unless they were specified
if isempty(opts.adducts)
    switch opts.polarity
        case 'Negative'
            opts.adducts = '"-H","+Cl"';
        case 'Positive'
            opts.adducts = '"+H","+Na","+K"';
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [json] = jsonParts(fn,opts)
% Contains the individual parts of the json file...

a = sprintf('{"MS_Analysis": {"Analyzer": "%s", ',opts.analyser);
b = sprintf('"Detector_Resolving_Power": {"Resolving_Power": %d, "mz": %d}, ',opts.rp,opts.rpMZ);
c = sprintf('"Ionisation_Source": "%s", ',opts.ionSource);
d = sprintf('"Polarity": "%s"}, ',opts.polarity);
e = sprintf('"Sample_Information": {"Condition": "%s", "Organism": "%s", ',opts.condition,opts.organism);
f = sprintf('"Organism_Part": "%s", ',opts.orgPart);
g = '"Sample_Growth_Conditions": "N/A"}, ';
h = '"Sample_Preparation": {"MALDI_Matrix": "N/A", ';
i = '"MALDI_Matrix_Application": "N/A", ';
j = sprintf('"Sample_Stabilisation": "%s", ',opts.sampStab);
k = '"Tissue_Modification": "N/A"}, ';
l = '"Submitted_By": {"Institution": "ICL", "Principal_Investigator": {"Email": "z.takats@imperial.ac.uk", "First_Name": "Zoltan", "Surname": "Takats"}, "Submitter": {"Email": "j.mckenzie@imperial.ac.uk", "First_Name": "James", "Surname": "McKenzie"}}, ';
m = sprintf('"metaspace_options": {"Dataset_Name": "%s", ',fn);
n = sprintf('"Metabolite_Database": [%s], ',opts.database);
o = sprintf('"Adducts": [%s]}, ',opts.adducts);
p = '"Additional_Information": {"Expected_Molecules_Freetext": "", ';
q = sprintf('"Sample_Description_Freetext": "%s"}}',opts.addInfo);

json = [a b c d e f g h i j k l m n o p q];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
