classdef emcStarFileReader
  %UNTITLED2 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = private)
    
    current_position_in_stack;
    current_column;

    position_in_stack_column;
    image_is_active_column;
    psi_column;
    theta_column;
    phi_column;
    x_shift_column;
    y_shift_column;
    defocus_1_column;
    defocus_2_column;
    defocus_angle_column;
    phase_shift_column;
    occupancy_column;
    logp_column;
    sigma_column;
    score_column;
    score_change_column;
    pixel_size_column;
    microscope_voltage_kv_column;
    microscope_spherical_aberration_mm_column;
    amplitude_contrast_column;
    beam_tilt_x_column;
    beam_tilt_y_column;
    image_shift_x_column;
    image_shift_y_column;
    stack_filename_column;
    original_image_filename_column;
    reference_3d_filename_column;
    best_2d_class_column;
    beam_tilt_group_column;
    particle_group_column;
    pre_exposure_column;
    total_exposure_column;
    
  end
  
  properties (Access = public)
  
    filename;
    input_file;

    using_external_array;

    % Pointer to an array of structs in cisTEM, use a cell
    cached_parameters;
    parameters_that_were_read;
    
  end
  
  
  
  methods
    
    function obj = emcStarFileReader(varargin)
      if isempty(varargin)
        % Do nothing
      else
        if isa(varargin{1}, 'string') && isa(alternate_cached_parameters_pointer, 'cell')
          obj.Init(varargin{1}, varargin{2}, varargin{3});
        else
          error('inputs to emcStarFileReader should be a filename, and cell');
        end
      end
    end
    
    function obj = Init(obj, wanted_filename, alternate_cached_parameters_pointer, exclude_negative_film_numbers)

    end

   
    function output =  ReturnPositionInStack(obj, line_number) ; output = obj.cached_parameters{line_number}.position_in_stack; end
    function output =  ReturnImageIsActive(obj, line_number) ; output = obj.cached_parameters{line_number}.image_is_active; end
    function output =  ReturnPhi(obj, line_number) ; output = obj.cached_parameters{line_number}.phi; end
    function output =  ReturnTheta(obj, line_number) ; output = obj.cached_parameters{line_number}.theta; end
    function output =  ReturnPsi(obj, line_number) ; output = obj.cached_parameters{line_number}.psi; end
    function output =  ReturnXShift(obj, line_number) ; output = obj.cached_parameters{line_number}.x_shift; end
    function output =  ReturnYShift(obj, line_number) ; output = obj.cached_parameters{line_number}.y_shift; end
    function output =  ReturnDefocus1(obj, line_number) ; output = obj.cached_parameters{line_number}.defocus_1; end
    function output =  ReturnDefocus2(obj, line_number) ; output = obj.cached_parameters{line_number}.defocus_2; end
    function output =  ReturnDefocusAngle(obj, line_number) ; output = obj.cached_parameters{line_number}.defocus_angle; end
    function output =  ReturnPhaseShift(obj, line_number) ; output = obj.cached_parameters{line_number}.phase_shift; end
    function output =  ReturnLogP(obj, line_number) ; output = obj.cached_parameters{line_number}.logp; end
    function output =  ReturnSigma(obj, line_number) ; output = obj.cached_parameters{line_number}.sigma; end
    function output =  ReturnScore(obj, line_number) ; output = obj.cached_parameters{line_number}.score; end
    function output =  ReturnScoreChange(obj, line_number) ; output = obj.cached_parameters{line_number}.score_change; end
    function output =  ReturnPixelSize(obj, line_number) ; output = obj.cached_parameters{line_number}.pixel_size; end
    function output =  ReturnMicroscopekV(obj, line_number) ; output = obj.cached_parameters{line_number}.microscope_voltage_kv; end
    function output =  ReturnMicroscopeCs(obj, line_number) ; output = obj.cached_parameters{line_number}.microscope_spherical_aberration_mm; end
    function output =  ReturnAmplitudeContrast(obj, line_number) ; output = obj.cached_parameters{line_number}.amplitude_contrast; end
    function output =  ReturnBeamTiltX(obj, line_number) ; output = obj.cached_parameters{line_number}.beam_tilt_x; end
    function output =  ReturnBeamTiltY(obj, line_number) ; output = obj.cached_parameters{line_number}.beam_tilt_y; end
    function output =  ReturnImageShiftX(obj, line_number) ; output = obj.cached_parameters{line_number}.image_shift_x; end
    function output =  ReturnImageShiftY(obj, line_number) ; output = obj.cached_parameters{line_number}.image_shift_y; end
    function output =  ReturnStackFilename(obj, line_number) ; output = obj.cached_parameters{line_number}.stack_filename; end
    function output =  ReturnOriginalImageFilename(obj, line_number) ; output = obj.cached_parameters{line_number}.original_image_filename; end
    function output =  ReturnReference3DFilename(obj, line_number) ; output = obj.cached_parameters{line_number}.reference_3d_filename; end
    function output =  ReturnBest2DClass(obj, line_number) ; output = obj.cached_parameters{line_number}.best_2d_class; end
    function output =  ReturnBeamTiltGroup(obj, line_number) ; output = obj.cached_parameters{line_number}.beam_tilt_group; end
    function output =  ReturnParticleGroup(obj, line_number) ; output = obj.cached_parameters{line_number}.particle_group; end
    function output =  ReturnPreExposure(obj, line_number) ; output = obj.cached_parameters{line_number}.pre_exposure; end
    function output =  ReturnTotalExpsosure(obj, line_number) ; output = obj.cached_parameters{line_number}.total_exposure; end
    
    [ ] = Open(obj,wanted_filename, alternate_cached_parameters_pointer);
    [ ] = ReadFile(obj,wanted_filename, error_string, alternate_cached_parameters_pointer, exclude_negative_film_number);
    [ ] = ExtractParametersFromLine(obj,wanted_line,error_string,  exclude_negative_film_numbers);

    [ ] = Close(obj);
    [ ] = Reset(obj);
    [ ] = ResetColumnPositions(obj);
   
  end
end

