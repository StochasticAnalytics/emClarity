function  EMC_assert_integer(input_val, varargin)

    % Use the default for any length
    assert_length = false;
    assert_passed = false;
    if ( nargin == 2 )
        assert_length = true;
        wanted_numel = varargin{1};
    elseif ( nargin > 2 )
        error('EMC_assert_numeric: too many input arguments');
    end

    if ( isa(input_val, 'integer') )
        if ( assert_length )
            if ( numel(input_val) == wanted_numel ) 
                assert_passed = true;
            end
        else
            assert_passed = true;
        end 
        return;
    end
    
    if ( assert_passed == false )
        error('EMC_assert_integer: input is not an integer or has wrong number of elements');
    end
    

end