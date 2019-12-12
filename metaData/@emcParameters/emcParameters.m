classdef emcParameters
  %Modelled after the cisTEMParameters class
  %   Detailed explanation goes here
  
  properties
    header_comments = "";
    all_parameters = {}; % array of parameter lines
    parameters_to_write; % binary mask 
    parameters_that_were_read; % binary mask
    
    % for defocus dependence
    average_defocus = 0.0; % Angstrom
    defocus_coeff_a;
    defocus_coeff_b;
  end
  
  methods
    
    % The constructure is default
    function obj = emcParameters()
      % Does nothing.
    end


    function output_long = ReturnNumberofLines(obj); output_long  = length(obj.all_parameters); end
    function output_emcParameterLine = ReturnLine(obj, line_number); output_emcParameterLine = obj.all_parameters{line_number}; end
    function output_int = ReturnPositionInStack(obj, line_number); output_int =  obj.all_parameters{line_number}.position_in_stack; end
    function output_int = ReturnImageIsActive(obj, line_number); output_int =  obj.all_parameters{line_number}.image_is_active; end
    function output_float =  ReturnPhi(obj, line_number); output_float =  obj.all_parameters{line_number}.phi; end
    function output_float =  ReturnTheta(obj, line_number); output_float =  obj.all_parameters{line_number}.theta; end
    function output_float =  ReturnPsi(obj, line_number); output_float =  obj.all_parameters{line_number}.psi; end
    function output_float =  ReturnXShift(obj, line_number); output_float =  obj.all_parameters{line_number}.x_shift; end
    function output_float =  ReturnYShift(obj, line_number); output_float =  obj.all_parameters{line_number}.y_shift; end
    function output_float =  ReturnDefocus1(obj, line_number); output_float =  obj.all_parameters{line_number}.defocus_1; end
    function output_float =  ReturnDefocus2(obj, line_number); output_float =  obj.all_parameters{line_number}.defocus_2; end
    function output_float =  ReturnDefocusAngle(obj, line_number); output_float =  obj.all_parameters{line_number}.defocus_angle; end
    function output_float =  ReturnPhaseShift(obj, line_number); output_float =  obj.all_parameters{line_number}.phase_shift; end
    function output_float =  ReturnOccupancy(obj, line_number); output_float =  obj.all_parameters{line_number}.occupancy; end
    function output_float =  ReturnLogP(obj, line_number); output_float =  obj.all_parameters{line_number}.logp; end
    function output_float =  ReturnSigma(obj, line_number); output_float =  obj.all_parameters{line_number}.sigma; end
    function output_float =  ReturnScore(obj, line_number); output_float =  obj.all_parameters{line_number}.score; end
    function output_float =  ReturnScoreChange(obj, line_number); output_float =  obj.all_parameters{line_number}.score_change; end
    function output_float =  ReturnPixelSize(obj, line_number); output_float =  obj.all_parameters{line_number}.pixel_size; end
    function output_float =  ReturnMicroscopekV(obj, line_number); output_float =  obj.all_parameters{line_number}.microscope_voltage_kv; end
    function output_float =  ReturnMicroscopeCs(obj, line_number); output_float =  obj.all_parameters{line_number}.microscope_spherical_aberration_mm; end
    function output_float =  ReturnAmplitudeContrast(obj, line_number); output_float =  obj.all_parameters{line_number}.amplitude_contrast; end
    function output_float =  ReturnBeamTiltX(obj, line_number); output_float =  obj.all_parameters{line_number}.beam_tilt_x; end
    function output_float =  ReturnBeamTiltY(obj, line_number); output_float =  obj.all_parameters{line_number}.beam_tilt_y; end
    function output_float =  ReturnImageShiftX(obj, line_number); output_float =  obj.all_parameters{line_number}.image_shift_x; end
    function output_float =  ReturnImageShiftY(obj, line_number); output_float =  obj.all_parameters{line_number}.image_shift_y; end
    function output_string = 	ReturnStackFilename(obj, line_number); output_string =  obj.all_parameters{line_number}.stack_filename; end
    function output_string =  ReturnOriginalImageFilename(obj, line_number); output_string =  obj.all_parameters{line_number}.original_image_filename; end
    function output_string =  ReturnReference3DFilename(obj, line_number); output_string =  obj.all_parameters{line_number}.reference_3d_filename; end
    function output_int =  ReturnBest2DClass(obj, line_number); output_int =  obj.all_parameters{line_number}.best_2d_class; end
    function output_int =  ReturnBeamTiltGroup(obj, line_number); output_int =  obj.all_parameters{line_number}.beam_tilt_group; end
    function output_int =  ReturnFrameNumber(obj, line_number); output_int =  obj.all_parameters{line_number}.particle_group; end
    function output_float =  ReturnPreExposure(obj, line_number); output_float =  obj.all_parameters{line_number}.pre_exposure; end
    function output_float =  ReturnTotalExposure(obj, line_number); output_float =  obj.all_parameters{line_number}.total_exposure; end

    [ output_float ] = ReturnAverageSigma(obj,exclude_negative_film_numbers);
    [ output_float ] = ReturnAverageOccupancy(obj,exclude_negative_film_numbers);
    [ output_float ] = ReturnAverageScore(obj,exclude_negative_film_numbers);

    [ ] = RemoveSigmaOutliers(obj,wanted_standard_deviation, exclude_negative_film_numbers, reciprocal_square);
    [ ] = RemoveScoreOutliers(obj,wanted_standard_deviation, exclude_negative_film_numbers, reciprocal_square);

    [ ] = CalculateDefocusDependence(obj,exclude_negative_film_numbers);
    [ ] = AdjustScores(obj,exclude_negative_film_numbers);
    [ output_float ] = ReturnScoreAdjustment(obj,defocus);
    [ output_float ] = ReturnScoreThreshold(obj,wanted_percentage, exclude_negative_film_numbers);

    [ output_float ] = ReturnMinScore(obj,exclude_negative_film_numbers);
    [ output_float ] = ReturnMaxScore(obj,exclude_negative_film_numbers);
    [ ouput_int ] = ReturnMinPositionInStack(obj,exclude_negative_film_numbers);
    [ ouput_int ] = ReturnMaxPositionInStack(obj,exclude_negative_film_numbers);

    [ ] = SetAllReference3DFilename(obj,wanted_filename);


    [ ] = ReadFromStarFile(obj,wanted_filename, program_name, exclude_negative_film_numbers );
    [ ] = AddCommentToHeader(obj,comment_to_add);
    [ ] = WriteTocisTEMStarFile(obj,wanted_filename, first_line_to_write, last_line_to_write, first_image_to_write, last_image_to_write);


    [ ] = PreallocateMemoryAndBlank(obj,number_to_allocate);
    [output_cisTEMParameterLine] = ReturnParameterVariances(obj,only_average_active);
    
    [ ] = ClearAll(obj)  
    [ ] = SortByReference3DFilename(obj);
  
  end
  




  
end

