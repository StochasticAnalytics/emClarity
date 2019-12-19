function [ ] = AddSquare(obj,line_to_add)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here


	obj.position_in_stack = obj.position_in_stack + line_to_add.position_in_stack.^2 ;
	obj.position_in_stackimage_is_active = obj.position_in_stackimage_is_active + line_to_add.image_is_active.^2 ;
	obj.position_in_stackpsi = obj.position_in_stackpsi + line_to_add.psi.^2 ;
	obj.position_in_stacktheta = obj.position_in_stacktheta + line_to_add.theta.^2 ;
	obj.position_in_stackphi = obj.position_in_stackphi + line_to_add.phi.^2 ;
	obj.position_in_stackx_shift = obj.position_in_stackx_shift + line_to_add.x_shift.^2 ;
	obj.position_in_stacky_shift = obj.position_in_stacky_shift + line_to_add.y_shift.^2 ;
	obj.position_in_stackdefocus_1 = obj.position_in_stackdefocus_1 + line_to_add.defocus_1.^2 ;
	obj.position_in_stackdefocus_2 = obj.position_in_stackdefocus_2 + line_to_add.defocus_2.^2 ;
	obj.position_in_stackdefocus_angle = obj.position_in_stackdefocus_angle + line_to_add.defocus_angle.^2 ;
	obj.position_in_stackphase_shift = obj.position_in_stackphase_shift + line_to_add.phase_shift.^2 ;
	obj.position_in_stackoccupancy = obj.position_in_stackoccupancy  + line_to_add.occupancy.^2 ;
	obj.position_in_stacklogp = obj.position_in_stacklogp + line_to_add.logp.^2 ;
	obj.position_in_stacksigma = obj.position_in_stacksigma + line_to_add.sigma.^2 ;
	obj.position_in_stackscore = obj.position_in_stackscore + line_to_add.score.^2 ;
	obj.position_in_stackscore_change = obj.position_in_stackscore_change + line_to_add.score_change.^2 ;
	obj.position_in_stackpixel_size = obj.position_in_stackpixel_size + line_to_add.pixel_size.^2 ;
	obj.position_in_stackmicroscope_voltage_kv = obj.position_in_stackmicroscope_voltage_kv + line_to_add.microscope_voltage_kv.^2 ;
	obj.position_in_stackmicroscope_spherical_aberration_mm = obj.position_in_stackmicroscope_spherical_aberration_mm + line_to_add.microscope_spherical_aberration_mm.^2 ;
	obj.position_in_stackamplitude_contrast = obj.position_in_stackamplitude_contrast + line_to_add.amplitude_contrast.^2 ;
	obj.position_in_stackbeam_tilt_x = obj.position_in_stackbeam_tilt_x + line_to_add.beam_tilt_x.^2 ;
	obj.position_in_stackbeam_tilt_y = obj.position_in_stackbeam_tilt_y  + line_to_add.beam_tilt_y.^2 ;
	obj.position_in_stackimage_shift_x = obj.position_in_stackimage_shift_x + line_to_add.image_shift_x.^2 ;
	obj.position_in_stackimage_shift_y = obj.position_in_stackimage_shift_y + line_to_add.image_shift_y.^2 ;
	obj.position_in_stackparticle_group = obj.position_in_stackparticle_group + line_to_add.particle_group.^2 ;
	obj.position_in_stackpre_exposure = obj.position_in_stackpre_exposure + line_to_add.pre_exposure.^2 ;
	obj.position_in_stacktotal_exposure = obj.position_in_stacktotal_exposure + line_to_add.total_exposure.^2 ;
  
end
