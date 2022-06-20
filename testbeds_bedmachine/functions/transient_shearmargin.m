function rheoB_weak = transient_shearmargin(rheoB, X, Y, left_mar, right_mar, amp, length)
%TRANSIENT_SHEARMARGIN This function creates a transient shear margin
%weakening. The weakenining is approximated by a parabolic function
%
%  Input:
%       rheoB        [double]: B number in the constitutive equation
%       X      [double array]: meshgrid X coordinates
%       Y      [double array]: meshgrid Y coordinates
%       left_mar  [int array]: vector of index representing left shear
%                              margin location
%       right_mar [int array]: vector of index representing right shear
%                              margin location
%       amp          [double]: amplitude, the lowerest point of parabola
%       length       [double]: length of the shear margin
%
    
    % a check first
    if numel(left_mar) ~= numel(right_mar)
        disp('Shear margin width not identical!')
        return
    end
    
    % width unit length must be an odd number
    if rem(numel(left_mar), 2) == 0
        disp('Shear margin width in units must be an odd number!')
        return
    end
    
    % first, make a map of constant rheoB
    rheoB_cons = rheoB*ones(size(X));
    rheoB_weak_min = amp*rheoB;
    % given grid is regular, find the across-flow spatial interval
    interval = Y(2,1) - Y(1,1);
    dist = (left_mar(end) - left_mar(1))*interval; % distance in meter
    % equation y = ax^2 + bx + c, where b = 0, c = rheoB_weak_min
    % equating y = rheoB_cons, the constant value, the distance between two
    % intersection should be 2*sqrt((d-c)/a), which should be cell dist
    c = rheoB_weak_min;
    d = rheoB;
    a = 4*(d-c)/dist^2;
    b = 0; % centered at zero. We make temp_y accordingly, symmetric around 0
    
    % make a temporary y coordinate where the lowest point is at zero
    s_margin_unitwid = numel(left_mar);
    one_side_unitwid = (s_margin_unitwid-1)/2;
    temp_y = transpose((-1)*interval*one_side_unitwid:interval:interval*one_side_unitwid);
    parabola = a.*temp_y.^2 + b.*temp_y + c;
    % repmat, the whole along flow direction
    parabola_all = repmat(parabola, 1, size(X, 2));
    
    % substitute in the original rheoB matrix
    rheoB_weak = rheoB_cons;
    rheoB_weak(left_mar,1:floor(length/150))  = parabola_all(:,1:floor(length/150));
    rheoB_weak(right_mar,1:floor(length/150)) = parabola_all(:,1:floor(length/150));

end

