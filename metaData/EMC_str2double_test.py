import numpy as np
import EMC_str2double as emc

def test_EMC_str2double():
    # Test valid input
    assert emc.EMC_str2double('3.14') == 3.14
    assert emc.EMC_str2double('-2.718') == -2.718
    assert emc.EMC_str2double('1e6') == 1e6
    assert emc.EMC_str2double('1.23e-4') == 1.23e-4
    assert emc.EMC_str2double('inf') == np.inf
    assert emc.EMC_str2double('-inf') == -np.inf
    assert np.isnan(emc.EMC_str2double('nan'))


    try:
        emc.EMC_str2double('not a number')
    except ValueError as e:
        assert str(e) == "EMC_str2double: input string is not a number!\nReceived: not a number"

# Run test
test_EMC_str2double()