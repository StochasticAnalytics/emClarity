import sys
from EMC_parseParameterFile import EMC_parseParameterFile

if __name__ == '__main__':
    # Check that the script was called with the correct number of arguments
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} PARAMETER_FILE")
        sys.exit(1)
    
    # Get the parameter file name from the command-line argument
    PARAMETER_FILE = sys.argv[1]
    
    # Call EMC_parseParameterFile with the parameter file name
    pStruct = EMC_parseParameterFile(PARAMETER_FILE)
    
    # Print the resulting parameter structure
    print(pStruct)