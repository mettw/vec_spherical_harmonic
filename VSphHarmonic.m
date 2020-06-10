classdef VSphHarmonic < handle
    % SphHarmonic - Vector Spherical Harmonics
    %
    % eg
    %
    % vharm = VSphHarmonic("dummy", "dummy", 10);
    % vharm.set_l_m_iseven(1,-1,Parity.Complex);
    % vharm.plot_MRadiation;
    properties
        % The default size of the vectors is too small to make sense of the
        % plots.
        arrow_scale = 6;
        
    end
    
    properties (SetAccess='private')
        % Group
        % Originally I wanted to initalise the object with the group and
        % irriducible representation and have VSphHarmonic choose the
        % correct values for l and m.  However, I am not sure that this is
        % so usefull anymore, so I will comment all of this code out.
        %group;
        % Irreducible Representation
        %irr_rep;
        
        %    parity   Whether or not the harmonic is even of odd.  That is,
        %              whether the scalar spherical harmonic has a cosine or
        %              sine dependence for phi.
        parity;
        %    l         degree of the spherical harmonic
        l;
        %    m         order of the spherical harmonic
        m;
        
        % How detailed the plot should be
        num_steps;
    end
    
    properties (GetAccess={?MSphHarmonic,?NSphHarmonic}, SetAccess='private', Hidden=true)
        %    Ymn       The scalar spherical harmonic defined over a cubical
        %              region of space.
        Yln;
        
        %    M spherical harmonic
        M;
        
        %    N spherical harmonic
        N;
        
        %    x,y,z     coordinates of the cubical region of space.
        x;
        y;
        z;
        
        % Cartesian coordianates of a sphere.
        x_sph;
        y_sph;
        z_sph;
    end
        
    properties (GetAccess='private', SetAccess='private', Hidden=true)
        %    Ymn       The scalar spherical harmonic defined over the 
        %              unit sphere
        Yln_sph;
        
        
        % Theta and Phi through the 3D volume of space
        Theta;
        Phi;
        R;
        
        % theta and phi values over a sphere.
        theta_sph;
        phi_sph;
        
        % Table of l, m and parity values for every group and irreducible
        % representation for 2D periodic structures.
        %l_m_iseven_values = ["dummy", "dummy", "0", "0", "even"];
    end
    
    methods
        %
        % SETUP functions
        %
        
        %function hObj = VSphHarmonic(group, irr_rep, num_steps)
        function hObj = VSphHarmonic(l, m, parity, num_steps)
            hObj.l = l;
            hObj.m = m;
            hObj.parity = parity;
            
            %hObj.group = group;
            %hObj.irr_rep = irr_rep;
            
            if nargin==3
                hObj.num_steps = 16;
            else
                hObj.num_steps = num_steps;
            end
            %{
            % Must be strings
            if ~isstring(group) && ischar(group)
                group = convertCharsToStrings(group);
            elseif ~isstring(group)
                ME = MException('VSphHarmonic:groupNotString', ...
                    'The value of group is not a string.');
                throw(ME)
            end
            if ~isstring(irr_rep) && ischar(irr_rep)
                group = convertCharsToStrings(irr_rep);
            elseif ~isstring(irr_rep)
                ME = MException('VSphHarmonic:irrrepNotString', ...
                    'The value of irr_rep is not a string.');
                throw(ME)
            end
            
            % For convenience we can look up the values of l and m and
            % whether or not the harmonic is odd or even from a table.
            %
            % To set the values manually use 'dummy' for the group and
            % irr_rep parameters and then call the set_l_m_iseven() method.
            
            % Find the group
            tmp = hObj.l_m_iseven_values(...
                hObj.l_m_iseven_values(:,1)==group, 2:5);
            if isempty(tmp)
                ME = MException('VSphHarmonic:illegalGroup', ...
                    'The group %s is not valid.', group);
                throw(ME)
            end
            % if l<0
            if str2num(tmp(1,2)) < 0 %#ok<ST2NM>
                ME = MException('VSphHarmonic:unsupportedGroup', ...
                    'The group %s is not supported yet.', group);
                throw(ME)
            end
            
            % Find the irreducible representation
            tmp = tmp(tmp(:,1)==irr_rep, 2:4);
            if isempty(tmp)
                ME = MException(['VSphHarmonic:'...
                    'illegalIrreducibleRepresentation'], ...
                    'The irreducible representation %s is not valid.', ...
                    irr_rep);
                throw(ME)
            end
            % if l<0
            if str2num(tmp(1,1)) < 0 %#ok<ST2NM>
                ME = MException(['VSphHarmonic:'...
                    'unsupportedIrreducibleRepresentation'], ...
                    ['The irreducible representation %s '...
                    'is not supported yet.'], irr_rep);
                throw(ME)
            end
            
            % Now assign the values
            hObj.l = str2num(tmp(1)); %#ok<ST2NM>
            hObj.m = str2num(tmp(2)); %#ok<ST2NM>
            if strcmp(tmp(3),'even')
                hObj.parity = Parity.Even;
            elseif strcmp(tmp(3),'odd')
                hObj.parity = Parity.Odd;
            else
                hObj.parity = Parity.Complex;
            end
            %}
            hObj.create_M_N();
        end
        %{
        % Set the values of these parameters directly.
        function set_l_m_iseven(hObj, l, m, parity)
            hObj.l = l;
            hObj.m = m;
            hObj.parity = parity;
            
            hObj.createYlm();
            hObj.M = MSphHarmonic(hObj);
            hObj.N = NSphHarmonic(hObj);
        end
        %}
        function set_num_steps(hObj, n)
            hObj.num_steps = n;
            hObj.create_M_N;
        end
        
        function set_arrow_scale(hObj, scale)
            hObj.arrow_scale = scale;
        end
        %% 
        %
        % PLOTTING functions
        %

        function plot_Ylm(hObj)
            % plot the spherical harmonic.
            Ylns = interp3(hObj.x,hObj.y,hObj.z,hObj.Yln,hObj.x_sph,hObj.y_sph,hObj.z_sph);
            %[Xr,Yr,Zr]=sph2cart(Phi2,Theta2-pi/2,real(Ymns).^2);
            figure;
            surf(hObj.x_sph.*real(Ylns).^2,hObj.y_sph.*real(Ylns).^2,hObj.z_sph.*real(Ylns).^2, real(Ylns),...
                    'EdgeColor', 'flat', 'FaceColor','interp');
            hObj.set_pbaspect();
            switch hObj.parity
                case Parity.Even
                    title_str = sprintf('Y_%d^{%d,c}', hObj.l, hObj.m);
                case Parity.Odd
                    title_str = sprintf('Y_%d^{%d,s}', hObj.l, hObj.m);
                otherwise
                    title_str = sprintf('Y_%d^%d', hObj.l, hObj.m);
            end
            %cb = colorbar;
            %ylabel(cb, ['Re\{' title_str '\}']);
            xlabel('x');
            ylabel('y');
            zlabel('z');
            title(title_str);
            colormap jet;
        end
        
        function plot_Mang(hObj)
            hObj.M.plot_ang();
        end
        
        function plot_MangVec(hObj)
            hObj.M.plot_ang_vec();
        end
        
        function plot_MRadiation(hObj)
            hObj.M.plot_ang_abs_vec();
        end
        
        function plot_Nang(hObj)
            hObj.N.plot_ang();
        end
        
        function plot_NangVec(hObj)
            hObj.N.plot_ang_vec();
        end
        
        function plot_NRadiation(hObj)
            hObj.N.plot_ang_abs_vec();
        end
    end
        %% 
        %
        % SUNDRY functions
        %
    methods (Access='private', Hidden=true)
        function create_M_N(hObj)
            % theta and phi values over a sphere.
            hObj.theta_sph = 0:pi/(2*hObj.num_steps):pi;
            hObj.phi_sph = 0:pi/hObj.num_steps:2*pi;

            % Cartesian coordianates of a sphere.
            hObj.x_sph = sin(hObj.theta_sph).*cos(hObj.phi_sph)';
            hObj.y_sph = sin(hObj.theta_sph).*sin(hObj.phi_sph)';
            hObj.z_sph = cos(hObj.theta_sph).*ones(size(hObj.phi_sph))';

            % Square region of space from -1 to 1 in all dimensions
            [hObj.x, hObj.y, hObj.z] = meshgrid(-1:2/hObj.num_steps:1,-1:2/hObj.num_steps:1,-1:2/hObj.num_steps:1);

            % Theta and Phi through the 3D volume of space
            [hObj.Phi, hObj.Theta, hObj.R] = cart2sph(hObj.x, hObj.y, hObj.z); 
            
            hObj.createYlm();
            hObj.M = MSphHarmonic(hObj);
            hObj.N = NSphHarmonic(hObj);
        end
        
        function createYlm(hObj)
            % Create the scalar spherical harmonic
            Lln=legendre(hObj.l,cos(hObj.Theta+pi/2));

            if hObj.l~=0
              Lln=squeeze(Lln(abs(hObj.m)+1,:,:,:));
            end

            a1=((2*hObj.l+1)/(4*pi));
            a2=factorial(hObj.l-abs(hObj.m))/factorial(hObj.l+abs(hObj.m));
            C=sqrt(a1*a2);
            switch hObj.parity
                case Parity.Even
                    hObj.Yln=C*Lln.*cos(hObj.m*hObj.Phi);
                case Parity.Odd
                    hObj.Yln=C*Lln.*sin(hObj.m*hObj.Phi);
                otherwise
                    hObj.Yln=C*Lln.*exp(-1i*hObj.m*hObj.Phi);
            end
            hObj.Yln_sph = interp3(hObj.x, hObj.y, hObj.z, hObj.Yln,...
                hObj.x_sph, hObj.y_sph, hObj.z_sph);
        end
        
        function set_pbaspect(hObj)
            % Sets the axes of the spherical harmonic plot so that each axis
            % has the same scale.  ie, none of the axes are cramped.
            switch hObj.l
                case 0
                    pba = [1 1 1];
                case 1
                    switch hObj.m
                        case -1
                            pba = [2 1 1];
                        case 0
                            pba = [1 1 2];
                        case 1
                            pba = [2 1 1];
                    end
                case 2
                    switch hObj.m
                        case -2
                            pba = [3 3 1];
                        case -1
                            pba = [3 1 3];
                        case 0
                            pba = [1 1 4];
                        case 1
                            pba = [3 1 3];
                        case 2
                            pba = [3 3 1];
                    end
                case 3
                    switch hObj.m
                        case -3
                            pba = [4 4 1];
                        case -2
                            pba = [3 3 2];
                        case -1
                            pba = [3 1 4];
                        case 0
                            pba = [1 1 6];
                        case 1
                            pba = [3 1 4];
                        case 2
                            pba = [3 3 2];
                        case 3
                            pba = [4 4 1];
                    end
                otherwise
                    pba = [1 1 1];
            end
            if hObj.l<=3
                pbaspect(pba);
            end
        end
    end
end

