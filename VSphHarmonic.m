classdef VSphHarmonic < handle
    % SphHarmonic - Vector Spherical Harmonics
    %
    % eg
    %
    % vharm = VSphHarmonic(1,-1,Parity.Complex, 10);
    % vharm.set_l_m_parity(1,-1,Parity.Complex);
    % vharm.plot_MRadiation;
    
    %properties (GetAccess='private', SetAccess='private')
    properties (SetAccess='private')
        %    parity   Whether or not the harmonic is even of odd.  That is,
        %              whether the scalar spherical harmonic has a cosine or
        %              sine dependence for phi.
        parity;
        %    l         degree of the spherical harmonic
        l;
        %    m         order of the spherical harmonic
        m;
        
        %    Ymn       The scalar spherical harmonic defined over a cubical
        %              region of space.
        Yln;
        % over the unit sphere
        Yln_sph;
        %    M spherical harmonic
        M;
        
        %    N spherical harmonic
        N;
        
        %    x,y,z     coordinates of the cubical region of space.
        x;
        y;
        z;
        num_steps;
        % Theta and Phi through the 3D volume of space
        Theta;
        Phi;
        R;
        
        % Cartesian coordianates of a sphere.
        x_sph;
        y_sph;
        z_sph;
        % theta and phi values over a sphere.
        theta_sph;
        phi_sph;
        
        % The default size of the vectors is too small to make sense of the
        % plots.
        arrow_scale = 6;
    end
    
    methods
        %
        % SETUP functions
        %
        
        function hObj = VSphHarmonic(l, m, parity, num_steps)
            
            if nargin==3
                hObj.num_steps = 16;
            else
                hObj.num_steps = num_steps;
            end
            
            if l < 0 || abs(m) > l
                ME = MException('VSphHarmonic: Illegal value for l or m');
                throw(ME)
            end
            
            
            % Now assign the values
            hObj.l = l;
            hObj.m = m;
            hObj.parity = parity;
            
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
        
        % Set the values of these parameters directly.
        function set_l_m_parity(hObj, l, m, parity)
            if l < 0 || abs(m) > l
                ME = MException('VSphHarmonic: Illegal value for l or m');
                throw(ME)
            end
            
            hObj.l = l;
            hObj.m = m;
            hObj.parity = parity;
            
            hObj.createYlm();
            hObj.M = MSphHarmonic(hObj);
            hObj.N = NSphHarmonic(hObj);
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
        %% 
        %
        % SUNDRY functions
        %

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

