function [] = SetAllToTrue(obj)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

	obj.position_in_stack = true;
	obj.image_is_active = true;
	obj.psi = true;
	obj.theta = true;
	obj.phi = true;
	obj.x_shift = true;
	obj.y_shift = true;
	obj.defocus_1 = true;
	obj.defocus_2 = true;
	obj.defocus_angle = true;
	obj.phase_shift = true;
	obj.occupancy = true;
	obj.logp = true;
	obj.sigma = true;
	obj.score = true;
	obj.score_change= true;
	obj.pixel_size = true;
	obj.microscope_voltage_kv = true;
	obj.microscope_spherical_aberration_mm = true;
	obj.amplitude_contrast = true;
	obj.beam_tilt_x = true;
	obj.beam_tilt_y = true;
	obj.image_shift_x = true;
	obj.image_shift_y = true;
	obj.stack_filename = true;
	obj.original_image_filename = true;
	obj.reference_3d_filename = true;
	obj.best_2d_class = true;
	obj.beam_tilt_group = true;
	obj.particle_group = true;
	obj.pre_exposure = true;
	obj.total_exposure = true;

end

