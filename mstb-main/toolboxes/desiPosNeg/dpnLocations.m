function [locn] = dpnLocations(layout)
%desiLocations - these are the locations for the various axes and labels
%and things in the desi window...


switch layout
    
    case 'simple'

        locn.opt.ax = [0.0333 0.51 0.45 0.45];        
        locn.ms1.ax = [0.0333 0.01 0.45 0.45];
        locn.ms2.ax = [0.5167 0.01 0.45 0.45];
        locn.fu.ax  = [0.5167 0.51 0.45 0.45];
        
        locn.sp.ax = [0.6 0.525 0.375 0.45];
        locn.mv.ax = [0.6 0.025 0.375 0.45];      
        locn.sp.vis = 'off';
        locn.mv.vis = 'off';

        locn.opt.lab1 = [0 0 0.01 0.01];
        locn.opt.lab2 = [0 0 0.01 0.01];

        locn.ms1.lab1 = [0 0 0.01 0.01];
        locn.ms1.lab2 = [0 0 0.01 0.01];
        
        locn.ms2.lab1 = [0 0 0.01 0.01];
        locn.ms2.lab2 = [0 0 0.01 0.01];
        
        
    case 'detail'
        
        locn.opt.ax = [0.025 0.51 0.3 0.45];
        locn.ms1.ax = [0.025 0.01 0.3 0.45];
        
        locn.ms2.ax = [0.35 0.01 0.3 0.45];        
        locn.fu.ax  = [0.35 0.51 0.3 0.45];
        
        locn.sp.ax  = [0.675 0.51 0.3 0.45];
        locn.mv.ax  = [0.675 0.01 0.3 0.45];
        
        locn.fu.vis = 'on';
        locn.sp.vis = 'on';
        locn.mv.vis = 'on';


        locn.opt.lab1 = [0 0 0.01 0.01];
        locn.opt.lab2 = [0 0 0.01 0.01];

        locn.ms1.lab1 = [0 0 0.01 0.01];
        locn.ms1.lab2 = [0 0 0.01 0.01];
        
        locn.ms2.lab1 = [0 0 0.01 0.01];
        locn.ms2.lab2 = [0 0 0.01 0.01];

        
        
end

end

