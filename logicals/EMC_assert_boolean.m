function  EMC_assert_boolean(input_val)

    if ( isa(input_val, 'logical') )
        return;
    else 
        if ( isa(input_val, 'numeric') )
            if ( (input_val == 0) || (input_val == 1) )
                return;
            else
                error('EMC_assert_boolean: input_val is numeric but not 0 or 1 to represent false or true');
            end
        end
    end

end