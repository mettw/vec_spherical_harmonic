classdef SumVSphHarmonic < handle
    % SphHarmonic - Vector Spherical Harmonics
    %
    % eg
    %
    % svh = SumVSphHarmonic([1 2], [0 2], [Parity.Even Parity.Odd], [true false], 16);
    % 
    
    %properties (GetAccess='private', SetAccess='private')
    properties (SetAccess='private')
        
        v_harmonics;
        Mangx;
        Mangy;
        Mangz;
        Mang;
        Nangx;
        Nangy;
        Nangz;
        Nang;
        MNangx;
        MNangy;
        MNangz;
        MNang;
        arrow_scale = 6;
        
    end
    
    methods
        %
        % SETUP functions
        %
        
        %function hObj = SumVSphHarmonic(group, irr_rep, num_steps)
        function hObj = SumVSphHarmonic(l_values, m_values, parity_vals, add_to_M, num_steps)
            hObj.Mangx = [];
            hObj.Mangy = [];
            hObj.Mangz = [];
            hObj.Mang = [];
            hObj.Nangx = [];
            hObj.Nangy = [];
            hObj.Nangz = [];
            hObj.Nang = [];
            
            if length(l_values) ~= length(m_values) || ...
                    length(l_values) ~= length(parity_vals) || ...
                    length(l_values) ~= length(add_to_M)
                ME = MException('SumVSphHarmonic:groupNotEqIrr_rep', ...
                    ['The vectors group, irr_rep and is_even and add_M must have the same',...
                    ' length [length(group) = %d, length(irr_rep) = %d,',...
                    ' length(is_even) = %d, length(add_M) = %d].'],...
                    length(l_values), length(m_values), length(parity_vals),...
                    length(add_to_M));
                throw(ME)
            end

            hObj.v_harmonics{1} = VSphHarmonic(l_values(1), m_values(1), parity_vals(1), num_steps);
            scale = 1;
            hObj.Mangx = scale*hObj.v_harmonics{1}.M.Mangx;
            hObj.Mangy = scale*hObj.v_harmonics{1}.M.Mangy;
            hObj.Mangz = scale*hObj.v_harmonics{1}.M.Mangz;
            hObj.Nangx = scale*hObj.v_harmonics{1}.N.Nangx;
            hObj.Nangy = scale*hObj.v_harmonics{1}.N.Nangy;
            hObj.Nangz = scale*hObj.v_harmonics{1}.N.Nangz;
            if add_to_M(1)
                hObj.MNangx = scale*hObj.v_harmonics{1}.M.Mangx;
                hObj.MNangy = scale*hObj.v_harmonics{1}.M.Mangy;
                hObj.MNangz = scale*hObj.v_harmonics{1}.M.Mangz;
            else
                hObj.MNangx = scale*hObj.v_harmonics{1}.N.Nangx;
                hObj.MNangy = scale*hObj.v_harmonics{1}.N.Nangy;
                hObj.MNangz = scale*hObj.v_harmonics{1}.N.Nangz;
            end
            for i=2:length(l_values)
                hObj.v_harmonics{i} = VSphHarmonic(l_values(i), m_values(i), parity_vals(i), num_steps);
                hObj.Mangx = hObj.Mangx + hObj.v_harmonics{i}.M.Mangx;
                hObj.Mangy = hObj.Mangy + hObj.v_harmonics{i}.M.Mangy;
                hObj.Mangz = hObj.Mangz + hObj.v_harmonics{i}.M.Mangz;
                hObj.Nangx = hObj.Nangx + hObj.v_harmonics{i}.N.Nangx;
                hObj.Nangy = hObj.Nangy + hObj.v_harmonics{i}.N.Nangy;
                hObj.Nangz = hObj.Nangz + hObj.v_harmonics{i}.N.Nangz;
                if add_to_M(i)
                    hObj.MNangx = hObj.MNangx + hObj.v_harmonics{i}.M.Mangx;
                    hObj.MNangy = hObj.MNangy + hObj.v_harmonics{i}.M.Mangy;
                    hObj.MNangz = hObj.MNangz + hObj.v_harmonics{i}.M.Mangz;
                else
                    hObj.MNangx = hObj.MNangx + hObj.v_harmonics{i}.N.Nangx;
                    hObj.MNangy = hObj.MNangy + hObj.v_harmonics{i}.N.Nangy;
                    hObj.MNangz = hObj.MNangz + hObj.v_harmonics{i}.N.Nangz;
                end
            end
            
            hObj.Mang = sqrt(hObj.Mangx.^2 + hObj.Mangy.^2 + hObj.Mangz.^2);
            hObj.Nang = sqrt(hObj.Nangx.^2 + hObj.Nangy.^2 + hObj.Nangz.^2);
            hObj.MNang = sqrt(hObj.MNangx.^2 + hObj.MNangy.^2 + hObj.MNangz.^2);
        end
        
        %% 
        %
        % PLOTTING functions
        %
        
        function plot_Mang(hObj)
            figure;

            surf(hObj.v_harmonics{1}.x_sph.*hObj.Mang, hObj.v_harmonics{1}.y_sph.*hObj.Mang, ...
                hObj.v_harmonics{1}.z_sph.*hObj.Mang, hObj.Mang,...
                'EdgeColor', 'flat', 'FaceColor','interp');
            title('$$|\mathrm{\sum (r/|r|)\times M}|$$','Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            colormap jet;
        end
        
        function plot_MangVec(hObj)
            figure;

            quiver3(hObj.v_harmonics{1}.x_sph, hObj.v_harmonics{1}.y_sph, hObj.v_harmonics{1}.z_sph, ...
                hObj.Mangx, hObj.Mangy, hObj.Mangz, hObj.v_harmonics{1}.arrow_scale,...
                'Color', 'black');
            pbaspect([1 1 1]);
            title('$$\mathrm{\sum (r/|r|)\times M}$$','Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            hold on;
            surf(hObj.v_harmonics{1}.x_sph, hObj.v_harmonics{1}.y_sph, ...
                hObj.v_harmonics{1}.z_sph, real(hObj.v_harmonics{1}.Yln_sph),...
                    'EdgeColor', 'interp', 'FaceColor','none');
            hold off;
            colormap jet;
        end
        
        function plot_MRadiation(hObj)
            figure;
            
            % The vector field is the real parts of M, so we want the
            % absolute value of the real parts of the vectors.
            abs_vec = sqrt(real(hObj.Mangx).^2 + real(hObj.Mangy).^2 + real(hObj.Mangz).^2);
            surf(hObj.v_harmonics{1}.x_sph.*abs_vec, hObj.v_harmonics{1}.y_sph.*abs_vec, hObj.v_harmonics{1}.z_sph.*abs_vec,...
                abs_vec, 'EdgeColor', 'flat', 'FaceColor','interp');
            
            pbaspect([1 1 1]);
            %hObj.set_pbaspect2;
            %{
            switch hObj.parent.parity 
                case Parity.Even
                    title_str = sprintf('$$\\mathrm{|(r/|r|)\\times M^e_{%d,%d}|}$$', hObj.parent.l, hObj.parent.m);
                case Parity.Odd
                    title_str = sprintf('$$\\mathrm{|(r/|r|)\\times M^o_{%d,%d}|}$$', hObj.parent.l, hObj.parent.m);
                otherwise
                    title_str = sprintf('$$\\mathrm{|(r/|r|)\\times M_{%d,%d}|}$$', hObj.parent.l, hObj.parent.m);
            end
            %}
            %title(title_str,'Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            colormap jet;
        end
        
        function plot_Nang(hObj)
            figure;

            surf(hObj.v_harmonics{1}.x_sph.*hObj.Nang, hObj.v_harmonics{1}.y_sph.*hObj.Nang, ...
                hObj.v_harmonics{1}.z_sph.*hObj.Nang, hObj.Nang,...
                'EdgeColor', 'flat', 'FaceColor','interp');
            title('$$|\mathrm{\sum (r/|r|)\times N}|$$','Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            colormap jet;
        end
        
        function plot_NangVec(hObj)
            figure;

            quiver3(hObj.v_harmonics{1}.x_sph, hObj.v_harmonics{1}.y_sph, hObj.v_harmonics{1}.z_sph, ...
                hObj.Nangx, hObj.Nangy, hObj.Nangz, hObj.v_harmonics{1}.arrow_scale,...
                'Color', 'black');
            pbaspect([1 1 1]);
            title('$$\mathrm{\sum (r/|r|)\times N}$$','Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            hold on;
            surf(hObj.v_harmonics{1}.x_sph, hObj.v_harmonics{1}.y_sph, ...
                hObj.v_harmonics{1}.z_sph, real(hObj.v_harmonics{1}.Yln_sph),...
                    'EdgeColor', 'interp', 'FaceColor','none');
            hold off;
            colormap jet;
        end
        
        function plot_MNang(hObj)
            figure;

            surf(hObj.v_harmonics{1}.x_sph.*abs(hObj.MNang), hObj.v_harmonics{1}.y_sph.*abs(hObj.MNang), ...
                hObj.v_harmonics{1}.z_sph.*abs(hObj.MNang), abs(hObj.MNang),...
                'EdgeColor', 'flat', 'FaceColor','interp');
            title('$$|\mathrm{\sum (r/|r|)\times N}|$$','Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            colormap jet;
        end
        
        function plot_MNangVec(hObj)
            figure;

            quiver3(hObj.v_harmonics{1}.x_sph, hObj.v_harmonics{1}.y_sph, hObj.v_harmonics{1}.z_sph, ...
                hObj.MNangx, hObj.MNangy, hObj.MNangz, hObj.v_harmonics{1}.arrow_scale,...
                'Color', 'black');
            pbaspect([1 1 1]);
            title('$$\mathrm{\sum (r/|r|)\times N}$$','Interpreter','latex');
            xlabel('x');
            ylabel('y');
            zlabel('z');
            hold on;
            surf(hObj.v_harmonics{1}.x_sph, hObj.v_harmonics{1}.y_sph, ...
                hObj.v_harmonics{1}.z_sph, real(hObj.v_harmonics{1}.Yln_sph),...
                    'EdgeColor', 'interp', 'FaceColor','none');
            hold off;
            colormap jet;
        end
        %% 
        %
        % SUNDRY functions
        %

    end
end

