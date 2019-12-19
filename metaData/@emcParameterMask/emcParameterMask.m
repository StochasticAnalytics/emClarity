classdef emcParameterMask
  %UNTITLED5 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    position_in_stack;
    image_is_active;
    psi;
    theta;
    phi;
    x_shift;
    y_shift;
    defocus_1;
    defocus_2;
    defocus_angle;
    phase_shift;
    occupancy;
    logp;
    sigma;
    score;
    score_change;
    pixel_size;
    microscope_voltage_kv;
    microscope_spherical_aberration_mm;
    amplitude_contrast;
    beam_tilt_x;
    beam_tilt_y;
    image_shift_x;
    image_shift_y;
    stack_filename;
    original_image_filename;
    reference_3d_filename;
    best_2d_class;
    beam_tilt_group;
    particle_group;
    pre_exposure;
    total_exposure;
    
  end
  
  methods
    function obj = emcParameterMask()
      % Do nothing
    end
    [output] = SetActiveParameters(obj,parameters_to_set); % uses takes the defines above bitwise
    SetAllToTrue(obj);
    SetAllToFalse(obj);
  end
end

