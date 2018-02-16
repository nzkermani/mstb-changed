function warp_pts = verb_plot_a2(msimage, warp_p, tmplt_pts, error_img, iterNum)
% VERB_PLOT_A - Verbose fitting plot
%   VERB_PLOT_A(V, WARP_P, TMPLT_PTS, ERROR_IMG)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: verb_plot_a.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin == 4
    iterNum = '0';
end

% Scaled error image
try
    set(msimage.op2msalign.herror,'CData',(error_img + 256) / 2);
    flag = true;
catch
    set(msimage(2),'CData',(error_img + 256) / 2);
    flag = false;
end

% title([num2str(min(100,5+ceil(100*sum((error_img(:)==0))./length(error_img(:))))),...
%     '% matched pixels.  (' int2str(iterNum) ')' char(10) ...
%     'Alignment in progress, please wait.'],...
%     'FontSize',10);

% Current parameters 
M = [warp_p; 0 0 1];
M(1,1) = M(1,1) + 1;
M(2,2) = M(2,2) + 1;
warp_pts =  M * [tmplt_pts; ones(1, size(tmplt_pts,2))];

% Draw the box in the appropriate axes
if flag
    set(msimage.op2msalign.haffine, ...
        'XData', [warp_pts(1,:) warp_pts(1,1)], ...
        'YData', [warp_pts(2,:) warp_pts(2,1)]); 
else
    set(msimage(1),...
        'XData', [warp_pts(1,:) warp_pts(1,1)], ...
        'YData', [warp_pts(2,:) warp_pts(2,1)]); 
end

% Flush refresh the drawing mechanism
drawnow;


end