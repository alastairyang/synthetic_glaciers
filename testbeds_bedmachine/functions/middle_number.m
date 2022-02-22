function mid_num = middle_number(num)
%MIDDLE_NUMBER Output the middle number
    
    if rem(num, 2) == 0 % it is even
        mid_num = num/2;
    else
        mid_num = (num+1)/2;
    end
    
end

