function [] = SetAllToFalse(obj)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

	obj.position_in_stack = false;
	obj.image_is_active = false;
	obj.psi = false;
	obj.theta = false;
	obj.phi = false;
	obj.x_shift = false;
	obj.y_shift = false;
	obj.defocus_1 = false;
	obj.defocus_2 = false;
	obj.defocus_angle = false;
	obj.phase_shift = false;
	obj.occupancy = false;
	obj.logp = false;
	obj.sigma = false;
	obj.score = false;
	obj.score_change= false;
	obj.pixel_size = false;
	obj.microscope_voltage_kv = false;
	obj.microscope_spherical_aberration_mm = false;
	obj.amplitude_contrast = false;
	obj.beam_tilt_x = false;
	obj.beam_tilt_y = false;
	obj.image_shift_x = false;
	obj.image_shift_y = false;
	obj.stack_filename = false;
	obj.original_image_filename = false;
	obj.reference_3d_filename = false;
	obj.best_2d_class = false;
	obj.beam_tilt_group = false;
	obj.particle_group = false;
	obj.pre_exposure = false;
	obj.total_exposure = false;
  
end

