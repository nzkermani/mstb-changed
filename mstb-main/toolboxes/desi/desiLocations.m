function [locn] = desiLocations(layout)
%desiLocations - these are the locations for the various axes and labels
%and things in the desi window...


switch layout
    
    case 'simple'

        locn.opt.ax = [0.025 0.2 0.45 0.7];
        locn.ms1.ax = [0.525 0.2 0.45 0.7];
        
        locn.sp.ax = [0.6 0.525 0.375 0.45];
        locn.mv.ax = [0.6 0.025 0.375 0.45];
        locn.sp.vis = 'off';
        locn.mv.vis = 'off';

        locn.opt.lab1 = [0.025 0.95 0.45 0.05];
        locn.opt.lab2 = [0.025 0.1 0.45 0.05];

        locn.ms1.lab1 = [0.525 0.95 0.45 0.05];
        locn.ms1.lab2 = [0.525 0.1 0.45 0.05];
        
        
    case 'detail'
        
        locn.opt.ax = [0.025 0.525 0.45 0.45];
        locn.ms1.ax = [0.025 0.025 0.45 0.45];
        
        locn.sp.ax = [0.525 0.525 0.45 0.45];
        locn.mv.ax = [0.525 0.025 0.45 0.45];
        locn.sp.vis = 'on';
        locn.mv.vis = 'on';


        locn.opt.lab1 = [0 0 0.01 0.01];
        locn.opt.lab2 = [0 0 0.01 0.01];

        locn.ms1.lab1 = [0 0 0.01 0.01];
        locn.ms1.lab2 = [0 0 0.01 0.01];

        
        
end

end

