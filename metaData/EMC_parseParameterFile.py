import EMC_str2double as emc

def EMC_parseParameterFile(PARAMETER_FILE):
    # Open the parameter file
    with open(PARAMETER_FILE, 'r') as fileID:
        # Read the file into a list of strings, ignoring comments and empty lines
        p = [line.strip() for line in fileID if line.strip() and not line.startswith('%')]
        
    # Check that all parameters are name: value pairs
    stringValues = ['subTomoMeta', 'Ali_mType', 'Cls_mType', 'Cls_mType', 'Raw_mType', 'Fsc_mType',
                    'Pca_distMeasure', 'Kms_mType', 'flgPrecision', 'Tmp_xcfScale', 'fastScratchDisk',
                    'Tmp_eraseMaskType', 'startingDirection', 'Peak_mType', 'symmetry']
    pStruct = {}
    for line in p:
        try:
            name, value = line.split('=', 1)
            name = name.strip()
            value = value.strip()
            # also strip any trailing ';' or ','
            if value[-1] in [';', ',']:
                value = value[:-1]
            if name in stringValues:
                pStruct[name] = value
            else:
                pStruct[name] = emc.EMC_str2double(value)
        except:
            err_msg = f"BH_parseParameterFile: invalid parameter line!\nReceived: {line}"
            raise ValueError(err_msg)
    
    return pStruct