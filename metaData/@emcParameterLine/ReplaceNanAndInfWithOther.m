function [ ] = ReplaceNanAndInfWithOther(obj,other_params)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

	if ~isfinite(obj.position_in_stack); obj.position_in_stack = other_params.position_in_stack; end
	if ~isfinite(obj.image_is_active); obj.image_is_active = other_params.image_is_active; end
	if ~isfinite(obj.psi); obj.psi = other_params.psi; end
	if ~isfinite(obj.theta); obj.theta = other_params.theta; end
	if ~isfinite(obj.phi); obj.phi = other_params.phi; end
	if ~isfinite(obj.x_shift); other_params.x_shift = other_params.x_shift; end
	if ~isfinite(obj.y_shift); obj.y_shift = other_params.y_shift; end
	if ~isfinite(obj.defocus_1); obj.defocus_1 = other_params.defocus_1; end
	if ~isfinite(obj.defocus_2); obj.defocus_2 = other_params.defocus_2; end
	if ~isfinite(obj.defocus_angle); obj.defocus_angle = other_params.defocus_angle; end
	if ~isfinite(obj.phase_shift); obj.phase_shift = other_params.phase_shift; end
	if ~isfinite(obj.occupancy); obj.occupancy = other_params.occupancy; end
	if ~isfinite(obj.logp); obj.logp = other_params.logp;end
	if ~isfinite(obj.sigma); obj.sigma = other_params.sigma; end
	if ~isfinite(obj.score); obj.score = other_params.score; end
	if ~isfinite(obj.score_change); obj.score_change = other_params.score_change; end
	if ~isfinite(obj.pixel_size); obj.pixel_size = other_params.pixel_size; end
	if ~isfinite(obj.microscope_voltage_kv); obj.microscope_voltage_kv = other_params.microscope_voltage_kv; end
	if ~isfinite(obj.microscope_spherical_aberration_mm); obj.microscope_spherical_aberration_mm = other_params.microscope_spherical_aberration_mm; end
	if ~isfinite(obj.amplitude_contrast); obj.amplitude_contrast = other_params.amplitude_contrast; end
	if ~isfinite(obj.beam_tilt_x); obj.beam_tilt_x  = other_params.beam_tilt_x ; end
	if ~isfinite(obj.beam_tilt_y); obj.beam_tilt_y = other_params.beam_tilt_y; end
	if ~isfinite(obj.image_shift_x); obj.image_shift_x = other_params.image_shift_x; end
	if ~isfinite(obj.image_shift_y); obj.image_shift_y = other_params.image_shift_y; end
	if ~isfinite(obj.stack_filename); obj.stack_filename = other_params.stack_filename; end
	if ~isfinite(obj.original_image_filename); obj.original_image_filename = other_params.original_image_filename ; end
	if ~isfinite(obj.reference_3d_filename); obj.reference_3d_filename = other_params.reference_3d_filename;
	if ~isfinite(obj.best_2d_class); obj.best_2d_class = other_params.best_2d_class;end
	if ~isfinite(obj.beam_tilt_group); obj.beam_tilt_group = other_params.beam_tilt_group;
	if ~isfinite(obj.particle_group); obj.particle_group = other_params.particle_group;end
	if ~isfinite(obj.pre_exposure); obj.pre_exposure = other_params.pre_exposure; end
	if ~isfinite(obj.total_exposure); obj.total_exposure  = other_params.total_exposure ;end
    
end
