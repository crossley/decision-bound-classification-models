function [a1,a2,b] = line_from_points(point1,point2)

% [a1,a2,b] = line_from_points(point1,point2)
% This function reads in two points and returns the coefficients 
% a1,a2 and b in the expression a1*x+a2*y+b = 0 that describes the line
% containing the points.  For the purpose of model fitting later we require
% that a1 is between -1 and 1 and a2 = sqrt(1-a1^2).

x1=point1(1);
y1=point1(2);
x2=point2(1);
y2=point2(2);

if (point1 == point2) 
    
    disp('The points are the same!')
    a1 = 0;
    a2 = 0;
    b =0;
    
elseif (x1 ~= x2)
    
    slope = (y2-y1)/(x2-x1);
    intercept = y1-slope*x1;

    a1 = -slope;
    a2 = 1;
    b = -intercept;
    
    const = sqrt(a1^2+a2^2);

    if (a2 < 0) 
        const = -const;
    end
           
    a1 = a1/const;
    a2 = a2/const;
    b = b/const;
    
else
   
    slope = (x2-x1)/(y2-y1);
    intercept = x1-slope*y1;

    a1 = 1;
    a2 = -slope;
    b = -intercept;
    
    const = sqrt(a1^2+a2^2);

    if (a2 < 0) 
        const = -const;
    end
           
    a1 = a1/const;
    a2 = a2/const;
    b = b/const;
    
end


       
