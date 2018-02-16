function [ res ] = compositionFinder( mass, numC, isoDist )
%compositionFinder This function accepts an accurate mass and an
%estimation of the number of carbons that is ultimately derived from the
%isotope cluster.  The procedure starts with the core elements (CHNOPS) and
%then expands to others, such as Cl,Br,I,Se,Na,K.

tStart = now;

% Need to define a distribution tolerance. This is the percentage over
% which it is acceptable for peaks to disagree.
ppmLevel = 10;

%% Start by defining accurate masses of the elements
table = [   12.0000000000;  1.00782503207;  6.015122795;    14.0030740048;  15.99491461956; 
            18.99840322;    22.9897692809;  23.985041700;   27.9769265325;  30.97376163; 
            31.97207100;    34.96885268;    38.96370668;    39.96259098;    74.9215965;
            73.9224764;     78.9183371;     126.904473];
elements = {'C'; 'H'; 'Li'; 'N'; 'O'; ...
    'F'; 'Na'; 'Mg'; 'Si'; 'P'; ...
    'S'; 'Cl'; 'K'; 'Ca'; 'As'; ...
    'Se'; 'Br'; 'I'};
numEl = numel(table);


%% Let's get some simple information from the isoDist
numC = round( (isoDist(2) / 1.1) * (100/isoDist(1)) );

% Scale the first peak to 100
%isoDist = 100 .* isoDist ./ isoDist(1);

%% So let us start to find some answers - we'll save the best few

% A storage structure for the results!
res = zeros(10, numEl+2);
res(:,numEl+1) = ppmLevel;

% Starting carbon count
nCarb = numC;%[17];

% Let's go!
for nc = nCarb
    
    nh = 0; nli = 0; nn = 0; no = 0;
    nf = 0; nna = 0; nmg = 0; nsi = 0; np = 0;
    ns = 0; ncl = 0; nk = 0; nca = 0; nas = 0;
    nse = 0; nbr = 0; ni = 0;    
    
    xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);    
    
    for ni = 0:ceil(xmass/table(18))
        xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for nbr = 0:ceil(xmass/table(17))
        xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for nse = 0:0%ceil(xmass/table(16))
        %xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for nas = 0:0%ceil(xmass/table(15))
        %xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for nca = 0:0%ceil(xmass/table(14))
        %xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for nk = 0:ceil(xmass/table(13))
        xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for ncl = 0:ceil(xmass/table(12))
        xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for ns = 0:ceil(xmass/table(11))
        xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for np = 0:ceil(xmass/table(10))
        xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);                    

    for nsi = 0:0%ceil(xmass/table(9))
        %xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for nmg = 0:0%ceil(xmass/table(8))
        %xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for nna = 0:ceil(xmass/table(7))
        xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for nf = 0:ceil(xmass/table(6))
        %xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for no = 0:ceil(xmass/table(5))
        xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for nn = 0:ceil(xmass/table(4))
        xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for nli = 0:0%ceil(xmass/table(3))
        %xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);

    for nh = 0:ceil(xmass/table(2))
        xmass = mass - ([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table);
                                                
        %% So we've reached the beginning of the actual processes involved in arriving at a composition!

        % calculate a metric to see if the proposal is sensible
        rdb = nc - (0.5*nh) + (0.5*nn) + 1;

        if rdb < 0 || rdb > nc
            % do nothing
        else

            % calculate the mass of this iteration
            zmass = [nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni] * table;

            % this is the mass difference
            dmass = 1e6 * (mass-zmass) / zmass;

            % is it better than the threshold?
            if abs(dmass) < ppmLevel

                % how about the calculation of the isotopic distribution!
                [iso] = isotopicCalculator([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni]);

                % let's compare the distribution to see if it is worth keeping!
                %minI = min([numel(iso) numel(isoDist)]);
                isoDiff = 0;                
                for r = 1:numel(isoDist)
                    try
                        isoDiff = isoDiff + (abs(iso(r) - isoDist(r)) / isoDist(r));
                    catch error                        
                        isoDiff = isoDiff + 1;
                    end
                end
                
                if isoDiff < 5 && abs(dmass) < (ppmLevel)
                    %figure;
                    %scatter(1:numel(iso), iso, 100, 'r', 'filled'); hold on;
                    %scatter(1:numel(isoDist), isoDist, 25, 'b', 'filled');

                    % so where is the worst one?
                    [~,rep] = max(abs(res(:,numEl+1))+res(:,numEl+2));

                    if numel(rep) > 1
                        error('More than one being replaced!');
                    end
                    rep = rep(1);

                    res(rep,1:numEl) = [nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni];
                    res(rep,numEl+1) = dmass;
                    res(rep,numEl+2) = isoDiff;

                    % display the composition of this iteration
                    %disp(['C' int2str(nc) ' H' int2str(nh) ' N' int2str(nn) ' O' int2str(no) ' P' int2str(np) ' S' int2str(ns) ' Cl' int2str(ncl) ' Br' int2str(nbr) char(9) 'm/z: ' num2str(zmass) char(9) 'ppm: (' int2str(round(dmass)) ')' char(9) int2str(rdb) char(9) 'isoDiff: ' num2str(isoDiff)]);
                    [string] = printCompound([nc nh nli nn no nf nna nmg nsi np ns ncl nk nca nas nse nbr ni], elements);
                    disp([string char(9) '%%%' char(9) 'm/z: ' num2str(zmass) char(9) 'ppm: (' int2str(round(dmass)) ')' char(9) int2str(rdb) char(9) 'isoDiff: ' num2str(isoDiff)]);
                    %disp([char(13)]);
                    %pause                                                                               
                end
            end

                                                                                                                                                    end
    end
    nh = 0;
    end
    nli = 0;
    end
    nn = 0;
    end
    no = 0;
    end
    nf = 0;
    end
    nna = 0;
    end
    nmg = 0;
    end
    nsi = 0;
    end
    np = 0;
    end
    ns = 0;
    end
    ncl = 0;
    end
    nk = 0;
    end
    nca = 0;
    end
    nas = 0;
    end
    nse = 0;
    end
    nbr = 0;
    end
    ni = 0;
end

%% Should do something fancy with the results!
figure;
plot(1:numEl, res(:,1:numEl));
set(gca, 'XTick', 1:numEl, 'XTickLabel', elements);

disp(['Time:' char(9) datestr(now-tStart, 'MM:SS.FFF')]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [string] = printCompound(comp, elements)
% A simple function to create a string that can be printed, and contains no
% mention of elements that don't appear whatsover...

numI = numel(comp);
string = ['--> ' char(9)];
for n = 1:numI
    
    if comp(n) > 0        
        string = [string elements{n} int2str(comp(n)) ' '];
    else
        % do nothing
    end
end

% No display the string!
%disp(string);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


