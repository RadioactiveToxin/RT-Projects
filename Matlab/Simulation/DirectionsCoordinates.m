classdef DirectionsCoordinates
    properties(Access= private)
        directions (1,:) {mustBeNumeric, mustBeFinite};
    end
    
    methods
        % Constructor
        function obj=DirectionsCoordinates(dir)
            if(nargin==1)
                obj.directions=dir;
            end
        end
        
        % Directions Probability: %
        function [East,West,North,South,Up,Down]=DirProb(obj)
            % e:=East w:=West n:=North s:=South
            e=obj.directions(1); w=obj.directions(2);
            n=obj.directions(3); s=obj.directions(4);
            East=num2cell([1:e]);
            West=num2cell([e+1:e+w]);
            North=num2cell([e+w+1:e+w+n]);
            South=num2cell([e+w+n+1:e+w+n+s]);
            % u:=Up d:=Down
            if(length(obj.directions)==6)
                u=obj.directions(5); d=obj.directions(6);
                Up=num2cell([e+w+n+s+1:e+w+n+s+u]);
                Down=num2cell([e+w+n+s+u+1:100]);
            end    
        end
        
        % Particle Coordinates:(x,y,z)
        function [x,y,z]=PartCoordi(obj,coordinates)
            % X,Y axes
            x=coordinates(1);
            y=coordinates(2);
            % Z axis
            if(length(coordinates)==3)
                z=coordinates(3);
            end
        end
        
    end
end