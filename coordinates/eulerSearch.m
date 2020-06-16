classdef eulerSearch < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    % Nothing to do until Init is called
    refine_top_N = 0;
    number_of_search_dimensions = 0;
    number_of_search_positions = 0;
    number_of_out_of_plane_angles = 0;
    number_of_angles_at_each_theta = [];
    best_parameters_to_keep = 0;
    list_of_search_parameters = {};      
    list_of_best_parameters = {};
    symmetry_symbol = 'C1';
    % Rotations are Z X Z', phi, theta, psi (azimuth, out of plane, in
    % plane)
    phi_max = 0.0; % 0 to 360
    phi_start = 0.0;
    theta_max = 0.0; % 0 to 180
    theta_step = 0.0;
    psi_max = 0.0; % 0 to 360
    psi_step = 0.0;
    resolution_limit = 0.0;
    max_search_x = 0.0;
    max_search_y = 0.0;
    bipolar_search = false;
    parameter_map = struct( ...
      'phi', [] , ...
      'theta', [], ...
      'psi', []);
    random_start_angle = false;
  end
  
  methods
    
    function [obj] = eulerSearch(wanted_symmetry_symbol, ...
                                 wanted_theta_max,...
                                 wanted_theta_step,...
                                 wanted_psi_max,...
                                 wanted_psi_step,...
                                 wanted_resolution_limit,...
                                 wanted_parameters_to_keep,...
                                 wanted_random_start_angle)
      
      if (wanted_theta_max < 0)
        wanted_theta_max = abs(wanted_theta_max);
        obj.bipolar_search = true;
      end
      
      obj.random_start_angle = wanted_random_start_angle;
                                
      obj.theta_max = wanted_theta_max;
      obj.theta_step = wanted_theta_step;
      % emClarity takes -angle:step:angle. This max is based on 0:360
      % though
      obj.psi_max = 2.*wanted_psi_max;
      obj.psi_step = wanted_psi_step;
      obj.symmetry_symbol = wanted_symmetry_symbol;

      SetSymmetryLimits(obj);
      CalculateGridSearchPositions(obj);

    end

    
  
    function [] = CalculateGridSearchPositions(obj)


      theta_max_local = obj.theta_max;
      obj.parameter_map.psi = -obj.psi_max./2 : obj.psi_step : obj.psi_max/2;
      obj.number_of_search_positions = 0; 
      
      theta_search =  [ 0 : obj.theta_step : theta_max_local ];
      if (obj.bipolar_search)
        theta_search = [theta_search, flip(180-theta_search)];
      end
      obj.parameter_map.theta = theta_search;
      obj.number_of_out_of_plane_angles = length(theta_search);
      obj.number_of_angles_at_each_theta = zeros(obj.number_of_out_of_plane_angles,1);
      obj.parameter_map.phi = cell(obj.number_of_out_of_plane_angles,1);
      

      obj.number_of_search_positions = 0;
      % Change this to include inplane angles explicitly.
      
      for iT = 1:obj.number_of_out_of_plane_angles
        theta = obj.parameter_map.theta(iT);
        
        if (theta == 0.0 || theta == 180.0)
          phi_step = obj.phi_max;
        else
          % angular sampling was adapted from Spider subroutine VOEA (Paul Penczek)
          phi_step = 1.*abs(obj.theta_step / sind(theta));
          if (phi_step > obj.phi_max) 
            phi_step = obj.phi_max;
          else
            phi_step = obj.phi_max / floor(obj.phi_max / phi_step + 0.5);
          end
        end
    
        if (obj.random_start_angle == true) 
          phi_start_local = phi_step / 2.0 * (rand(1) - 0.5);
        else
          phi_start_local = 0.0;
        end
	             
        obj.parameter_map.phi{iT} = [0:phi_step:obj.phi_max - 1] + phi_start_local;        

        obj.number_of_angles_at_each_theta(iT) = length(obj.parameter_map.phi{iT}) .* length(obj.parameter_map.psi);
      end
      
      
    end

    function [] = SetSymmetryLimits(obj)
    
      switch obj.symmetry_symbol(1)
         case 'C'
           if (length(obj.symmetry_symbol) < 2)
             error('Cyclic symmetry requires an int specifying CX');
           end

          obj.psi_max = min(obj.psi_max,360.0 / str2double(obj.symmetry_symbol(2:end)));
          obj.theta_max = min(180.0, obj.theta_max);
          obj.phi_max = 360.0; % This will be incompatible with "symmetry_constrained_search" in BH_mutli_gridAngleSEarch (or whatever)
          
        case 'D'
          % FIXME is this right?
        if ((length(obj.symmetry_symbol) < 2))
          error('D symmetry requires an int specifying DX');
        end
          obj.psi_max = min(360.0 / str2double(obj.symmetry_symbol(2:end)));
          obj.theta_max = min(obj.theta_max,90.0);
          obj.phi_max = 360.0;

         case 'O'
           if ((length(obj.symmetry_symbol) > 1))
             error('Octahedral symmetry requires no int');
           end
           
          obj.psi_max = min(obj.psi_max,90.0);
          obj.theta_max = min(obj.theta_max,90); % cisTEM uses 54.7 but I'm not sure why. (Oddly, This is also the magic angle for SS-NMR?)
          obj.phi_max = 90.0;
           
         case 'I'
           % Double check convention: TODO
           % 2 fold on Z, 5 fold 31.17 deg around X on Y axis, 3 fold 20.91
           % deg around Y on X. For I2 the X/Y axes are flipped
          if ((length(obj.symmetry_symbol) < 2))
            obj.psi_max = min(obj.psi_max,180.0);
            obj.theta_max = 31.7;
            obj.phi_max = 180.0;
          elseif strcmp(obj.symmetry_symbol,'2')
            obj.psi_max = min(obj.psi_max,180.0);
            obj.theta_max = 31.7;
            obj.phi_max = 180.0;           
          else
            error('Icosohedral can be I or I2, not (%s)',obj.symmetry_symbol);
          end
        otherwise
          error('symmetry symbol (%s) not recognized', obj.symmetry_symbol);
      end
    end


  end
end

