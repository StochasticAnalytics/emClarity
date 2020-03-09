classdef EMC_Fourier < handle
% Fourier transformer for EMC functions.
% This class supports half (non-redundant) and
% full (redundant) Discrete Fourier Transforms (DFT).
%

% Created: 8Mar2020, R2019a
% Version: v.1.0
%

    % These attributes cannot be changed from outside the class methods.
    properties (SetAccess = private)
        size_real 	% size of the real space image
        size_freq  	% size of the frequency space spectrum.
        method    	% 'gpu' or 'cpu'
        precision 	% 'single' or 'double'
        is3d
        isEven
        isOnGpu
        centered  	% change frequency order of the spectrums, bandpass and indexes.
        half      	% compute half transforms, bandpass and indexes.
        half_wrap

        % These are used by the getters to know if the indexes are already set.
        index_fftshift_isSet = false;
      	index_ifftshift_isSet = false;
        index_half2full_isSet = false;
    end

    % These attributes cannot be changed from outside the class methods.
    % These attributes depends on the others attributes defined above.
    properties (Dependent = true, SetAccess = private)
        bandpass
        index_fftshift
        index_ifftshift
        index_half2full
    end

    % Constructor
    methods
        function obj = EMC_Fourier(SIZE, METHOD, OPTION)
        %
        % obj = EMC_Fourier(SIZE, METHOD, OPTION)
        % Constructor of EMC_Fourier.
        %
        % Input:
        %   SIZE (vector):              Size (in pixel) of the 3d/2d grids to transform; [x, y, z] or [x, y].
        %                               NOTE: [1,1], [N,1] or [1,N] are not allowed.
        %                               NOTE: This cannot be changed after initialization.
        %
        %	METHOD (str):               Method to use; 'gpu' or 'cpu'.
        %                               NOTE: This can be changed after initialization using to().
        %
        %   OPTION (cell|struct):       Optional parameters.
        %                               If cell: {field,value ; ...}, note the ';' between parameters.
        %                               NOTE: Can be empty.
        %                               NOTE: Unknown fields will raise an error.
        %
        %	  -> 'precision' (str):     Precision of the transforms; 'single' or 'double'.
        %                               NOTE: This can be changed after initialization using to().
        %                               default = 'single'
        %
        %     -> 'half' (bool):         Wheter or not the transformer should compute (and expect) the
        %                               non-redundant spectrum/bandpass. 
        %                               NOTE: This cannot be changed after initialization.
        %                               default = true
        %
        %     -> 'centered' (bool): 	Wheter or not the spectrums/bandpass should be centered (fftshift).
        %                               This will automatically fftshift the transforms after fft and
        %                               ifftshift the transforms before ifft.
        %                               NOTE: This cannot be changed after initialization.
        %                               default = false
        %
        %     -> 'plans' (str):         'estimate', 'measure', 'patient', 'exhaustive', 'hybrid' or 'none'
        %                               Create the fftw plans; equivalent to fftw('planned', plans).
        %                               default = 'none'
        %
        % Output:
        %   obj (handle):               Instance of the EMC_Fourier class.
        %
        % Other EMC-files required:
        %   EMC_is3d, EMC_getOption
        %

        [obj.is3d, obj.size_real] = EMC_is3d(SIZE);
        obj.isEven = ~mod(obj.size_real, 2);

        % method
        if strcmpi(METHOD, 'gpu')
            obj.method = 'gpu';
            obj.isOnGpu = true;
        elseif strcmpi(METHOD, 'cpu')
            obj.method = 'cpu';
            obj.isOnGpu = false;
        else
            error('EMC:Fourier', "METHOD should be 'gpu' or 'cpu'")
        end

        OPTION = EMC_getOption(OPTION, {'precision', 'mode_half', 'mode_fftshift', 'plans'}, false);

        % precision
        if isfield(OPTION, 'precision')
            if ~(ischar(OPTION.precision) || isstring(OPTION.precision)) || ...
               ~strcmpi(OPTION.precision, 'single') && ~strcmpi(OPTION.precision, 'double')
                error('EMC:Fourier', "OPTION.precision should be 'single' or 'double'")
            end
            obj.precision = OPTION.precision;
        else
            obj.precision = 'single';  % default
        end

        % centered
        if isfield(OPTION, 'centered')
            if ~isscalar(OPTION.centered) && ~islogical(OPTION.centered)
                error('EMC:Fourier', 'OPTION.centered should be true or false')
            end
            obj.centered = OPTION.centered;
        else
            obj.centered = true;  % default
        end

        % half
        if isfield(OPTION, 'half')
            if ~isscalar(OPTION.half) && ~islogical(OPTION.half)
                error('EMC:Fourier', 'OPTION.half should be true or false')
            end
            obj.half = OPTION.half;
            if obj.half
                obj.size_freq = [floor(SIZE(1) / 2) + 1, SIZE(2:end)];
                if obj.centered; obj.half_wrap = 'c2nc'; else; obj.half_wrap = 'nc2nc'; end
            else
                obj.size_freq = obj.size_real;
            end
        else
            obj.half = true;  % default
            obj.size_freq = [floor(SIZE(1) / 2) + 1, SIZE(2:end)];
            obj.half_wrap = 'nc2nc';
        end

        % plans
        if ~isfield(OPTION, 'plans') || strcmpi(OPTION.plans, 'none')
            OPTION.plans = 'none';  % default
        else
            fftw('planned', OPTION.plans)
            fftn(zeros(obj.size_real));
        end

        end  % EMC_Fourier
    end

    methods (Access = public)
        DFT = fft(obj, IMAGE);
        IMAGE = ifft(obj, DFT);
        obj = setBandpass(obj, ARRAY, PIXEL, HIGHPASS, LOWPASS, OPTION);
        ARRAY = applyBandpass(obj, ARRAY, OPTION);
        obj = to(obj, STR);
    end
 	methods (Static)
        INDEX = getIndex(TYPE, SIZE, METHOD, OPTION);
        FULL = half2full(WRAP, HALF, SIZE);
    end

    % Setters
    methods
        function value = get.index_fftshift(obj)
            if ~obj.index_fftshift_isSet
                obj.index_fftshift_isSet = true;
                obj.index_fftshift = obj.getIndex('fftshift', obj.size_real, obj.method, ...
                                                  {'half', obj.half});
            end
            value = obj.index_fftshift;
        end
        function value = get.index_ifftshift(obj)
            if ~obj.index_ifftshift_isSet
                obj.index_ifftshift_isSet = true;
                obj.index_ifftshift = obj.getIndex('ifftshift', obj.size_real, obj.method, ...
                                                   {'half', obj.half});
            end
            value = obj.index_ifftshift;
        end
        function value = get.index_half2full(obj)
            if ~obj.index_half2full_isSet
                obj.index_half2full_isSet = true;
                obj.index_half2full = obj.getIndex(obj.half_wrap, obj.size_real, obj.method, ...
                                                  {'half', obj.half});
            end
            value = obj.index_half2full;
        end
    end
end
