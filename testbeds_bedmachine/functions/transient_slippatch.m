function fric_coef_slip = transient_slippatch(fric_coef, X, Y, loc, width, amp)
%TRANSIENT_SHEARMARGIN This function creates a transient 2-D gaussian
%slippery patch with given location in the bed, width, and amplitude
%
%  Input:
%       fric_coef    [double]: uniform frictional coefficient
%       X      [double array]: meshgrid X coordinates
%       Y      [double array]: meshgrid Y coordinates
%       loc    [double array]: indices of gaussian center x,y coordinate in X,Y
%       width        [double]: meter, peak width of the gaussian
%       amp          [double]: meter, amplitude (negative for slippery
%                              patch)

    
    % first, make a map of constant frictional coefficient
    fric_coef_cons = fric_coef*ones(size(X));
    % given grid is regular, find the across-flow spatial interval
    x0_i = loc(1);
    y0_i = loc(2);
    x0 = X(1, x0_i);
    y0 = Y(y0_i, 1);
    gaussian = amp.*exp(-((X-x0).^2/(2*width^2)+...
                          (Y-y0).^2/(2*width^2)));
    
    % add together
    fric_coef_slip = fric_coef_cons + gaussian;

end

