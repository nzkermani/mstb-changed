function pdPolaritySwitch(~,~,fig)
%pdPolaritySwitch - change ions between pos and neg modes

% Get the current text string of the button
cur = fig.polarity.String;

switch cur
    case 'NEG'
        % Need to switch to pos
        ions = ionsPOS;
        fig.polarity.String = 'POS';
        
    case 'POS'
        % Need to switch to neg
        ions = ionsNEG;
        
        fig.polarity.String = 'NEG';
        
end

fig.ionList.String = ions;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [i] = ionsNEG

i = {'255.2329';'281.2485';'465.3043';'536.5047'; '744.5548'; '885.5498'};

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [i] = ionsPOS

i = ['2345' char(10) 'sdfghn' char(10) 'sdfgh'];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%