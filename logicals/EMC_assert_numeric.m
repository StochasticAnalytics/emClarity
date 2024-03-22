function  EMC_assert_numeric(input_val, varargin)

    % Use the default for any length
    assert_length = false;
    assert_passed = false;
    assert_range = false;
    
    if ( nargin > 1 )
        assert_length = true;
        wanted_numel = varargin{1};
        if ~isa(wanted_numel, "numeric") 
            error('EMC_assert_numeric: second argument must be an integer');
        end
        if (nargin > 2)
            assert_range = true;
            range = varargin{2};
            if ( ~isa(range, 'numeric') || numel(range) ~= 2 )
                error('EMC_assert_numeric: third argument must be numeric with two elements');
            end
        end
        if ( nargin > 3 ) 
            error('EMC_assert_numeric: too many input arguments');
        end
    end

    % If it is not numeric, throw an error
    if ( isa(input_val, 'numeric') )
        % Does it have enough values?
        if ( assert_length )
            if ( numel(input_val) == wanted_numel ) 
                assert_passed = true;
            end
            if (assert_range)
                for val_in_range = 1:wanted_numel
                    if ( input_val(val_in_range) < range(1) || input_val(val_in_range) > range(2) )
                        assert_passed = false;
                    end
                end
            end
        else
            assert_passed = true;
        end 
    end
    
    if ( assert_passed == false )
        error('EMC_assert_numeric: input is not numeric or has wrong number of elements or is out of range');
    end
    

end