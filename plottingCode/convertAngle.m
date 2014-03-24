function [theta] = convertAngle(r1,r2,angle)
    if (angle == pi)
      theta = pi;
    else
      theta = atan(tan(angle)*(r1/r2));
      if (angle > pi/2 && angle <= 3*pi/2)
        theta = theta + pi;
      else if (angle > 3*pi/2)
        theta = theta + 2*pi;
      end
    end
end
