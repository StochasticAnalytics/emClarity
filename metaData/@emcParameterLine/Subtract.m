function [ ] = Subtract(obj,line_to_subtract)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

	obj.position_in_stack = obj.position_in_stack - line_to_subtract.position_in_stack;
	obj.position_in_stackimage_is_active = obj.position_in_stackimage_is_active - line_to_subtract.image_is_active;
	obj.position_in_stackpsi = obj.position_in_stackpsi - line_to_subtract.psi;
	obj.position_in_stacktheta = obj.position_in_stacktheta - line_to_subtract.theta;
	obj.position_in_stackphi = obj.position_in_stackphi - line_to_subtract.phi;
	obj.position_in_stackx_shift = obj.position_in_stackx_shift - line_to_subtract.x_shift;
	obj.position_in_stacky_shift = obj.position_in_stacky_shift - line_to_subtract.y_shift;
	obj.position_in_stackdefocus_1 = obj.position_in_stackdefocus_1 - line_to_subtract.defocus_1;
	obj.position_in_stackdefocus_2 = obj.position_in_stackdefocus_2 - line_to_subtract.defocus_2;
	obj.position_in_stackdefocus_angle = obj.position_in_stackdefocus_angle - line_to_subtract.defocus_angle;
	obj.position_in_stackphase_shift = obj.position_in_stackphase_shift - line_to_subtract.phase_shift;
	obj.position_in_stackoccupancy = obj.position_in_stackoccupancy  - line_to_subtract.occupancy;
	obj.position_in_stacklogp = obj.position_in_stacklogp - line_to_subtract.logp;
	obj.position_in_stacksigma = obj.position_in_stacksigma - line_to_subtract.sigma;
	obj.position_in_stackscore = obj.position_in_stackscore - line_to_subtract.score;
	obj.position_in_stackscore_change = obj.position_in_stackscore_change - line_to_subtract.score_change;
	obj.position_in_stackpixel_size = obj.position_in_stackpixel_size - line_to_subtract.pixel_size;
	obj.position_in_stackmicroscope_voltage_kv = obj.position_in_stackmicroscope_voltage_kv - line_to_subtract.microscope_voltage_kv;
	obj.position_in_stackmicroscope_spherical_aberration_mm = obj.position_in_stackmicroscope_spherical_aberration_mm - line_to_subtract.microscope_spherical_aberration_mm;
	obj.position_in_stackamplitude_contrast = obj.position_in_stackamplitude_contrast - line_to_subtract.amplitude_contrast;
	obj.position_in_stackbeam_tilt_x = obj.position_in_stackbeam_tilt_x - line_to_subtract.beam_tilt_x;
	obj.position_in_stackbeam_tilt_y = obj.position_in_stackbeam_tilt_y  - line_to_subtract.beam_tilt_y;
	obj.position_in_stackimage_shift_x = obj.position_in_stackimage_shift_x - line_to_subtract.image_shift_x;
	obj.position_in_stackimage_shift_y = obj.position_in_stackimage_shift_y - line_to_subtract.image_shift_y;
	obj.position_in_stackparticle_group = obj.position_in_stackparticle_group - line_to_subtract.particle_group;
	obj.position_in_stackpre_exposure = obj.position_in_stackpre_exposure - line_to_subtract.pre_exposure;
	obj.position_in_stacktotal_exposure = obj.position_in_stacktotal_exposure - line_to_subtract.total_exposure;
end
