from numpy import ones

def EMC_str2double(input_str):
    try:
        output_double = float(input_str)
    except ValueError:
        try:
            output_double = eval(input_str)
            if output_double is None:
                raise ValueError
        except ValueError:
            try:
                # check for ones() is in the string and use the numpy function to replace it
                if 'ones(' in input_str:
                    input_str = input_str.replace('ones(', 'ones((')
                    input_str = input_str.replace(')', ', dtype=float)')
                    print(input_str)
            except ValueError:
                err_msg = f"EMC_str2double: input string is not a number!\nReceived: {input_str}"
                raise ValueError(err_msg)
    return output_double