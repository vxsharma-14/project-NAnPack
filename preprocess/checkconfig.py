#**************************************************************************
def CheckSections(config, InFileName):
    '''Check that all sections exist in the configuration file.

    If not, raise an exception.
    '''
    sections = ['SETUP', 'DOMAIN', 'MESH', 'IC', 'BC', 'CONST', 'STOP',\
                'OUTPUT']
    for section in sections:
        if config.has_section(section) == True:
            print(f'Checking section {section}: Completed.')
        else:
            raise Exception(f'Section "{section}" not found.\nin file :\
 "{InFileName}".\nCheck configuration file.')

#**************************************************************************
def CheckBCSections(config, InFileName):
    '''Check that all sections exist in the configuration file.

    If not, raise an exception.
    '''
    sections = ['INLET', 'WALL', 'FAR-FIELD', 'OUTLET']
    for section in sections:
        if config.has_section(section) == True:
            print(f'Checking section {section}: Completed.')
        else:
            raise Exception(f'Section "{section}" not found.\nin file :\
 "{InFileName}".\nCheck BC configuration file.')

#**************************************************************************
def CheckSetupSection(config, State, Model, Dimension, InFileName):
    '''Check if all the options in section SETUP exist and the values are
    entered only from the available choices.

    If not, raise an exception.
    '''
    import preprocess.fetchoptions as opt
    fetch = opt.FetchOptions()
    
    state_options = fetch.StateOptions()
    model_options = fetch.ModelOptions()
    dim_options = fetch.DimensionOptions()
    

    if not State.upper() in state_options:
        raise Exception(f'Invalid STATE info: {State}.\nin file:\
 "{InFileName}".\nin section: SETUP.\nAvailable options are:\
 {state_options}.')

    if not Model.upper() in model_options:
        raise Exception(f'Invalid MODEL info: {Model}.\nin file:\
 "{InFileName}".\nin section: SETUP.\nAvailable options are:\
 {model_options}.')

    if not Dimension.upper() in dim_options:
        raise Exception(f'Invalid DIMENSION info: {Dimension}.\nin file:\
 "{InFileName}".\nin section: SETUP.\nAvailable options are:\
 {dim_options}.')

    print("User inputs in SETUP section check: Completed.")
#**************************************************************************
