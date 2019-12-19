function [ ] = SetAllToZero(obj)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here


	obj.position_in_stack = 0;
	obj.image_is_active = 0;
	obj.psi = 0.0;
	obj.theta = 0.0;
	obj.phi = 0.0;
	obj.x_shift = 0.0;
	obj.y_shift = 0.0;
	obj.defocus_1 = 0.0;
	obj.defocus_2 = 0.0;
	obj.defocus_angle = 0.0;
	obj.phase_shift = 0.0;
	obj.occupancy = 0.0;
	obj.logp = 0.0;
	obj.sigma = 0.0;
	obj.score = 0.0;
	obj.score_change = 0.0;
	obj.pixel_size = 0.0;
	obj.microscope_voltage_kv = 0.0;
	obj.microscope_spherical_aberration_mm = 0.0;
	obj.amplitude_contrast = 0.0;
	obj.beam_tilt_x = 0.0;
	obj.beam_tilt_y = 0.0;
	obj.image_shift_x = 0.0;
	obj.image_shift_y = 0.0;
	obj.stack_filename = ' ';
	obj.original_image_filename = ' ';
	obj.reference_3d_filename = ' ';
	obj.best_2d_class = 0;
	obj.beam_tilt_group = 0;
	obj.particle_group = 0;
	obj.pre_exposure = 0.0;
	obj.total_exposure = 0.0;
end
