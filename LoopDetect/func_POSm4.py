#function definitions for examples for LoopDetect

def func_POSm4(x,klin,knonlin):
    '''
	example function: chain model with positive feedback regulation as from Baum et al., 2016

    ::Parameters:

        - x (numpy array, list, tuple): numerical, ordered values of length 4; these are the 
          variable values.
        - klin: (numpy array, list, tuple): numerical, ordered values of length 8; these are 
          some of the kinetic parameters of the model
        - knonlin: 

    :Returns:

        A numpy array of length 4 that contains the time-derivatives of the four variables.

    :Details:

        It is a function delivering an ordinary differential equation system with chain structure
        with 4 species and a positive feedback regulation from the last species on the
        conversion between species 1 and 2. The system can give rise to oscillations; 
        see Baum et al., 2016 for details.

    :Reference:
    
        Baum K, Politi AZ, Kofahl B, Steuer R, Wolf J: Feedback, Mass Conservation and Reaction 
        Kinetics Impact the Robustness of Cellular Oscillationss. Plos Comput Biol. 
        2016;12(12):e1005298.
    '''
    dx = np.zeros(4)
    dx[0] = klin[0]-(klin[1]*(1 + x[3]/pow(knonlin[0],knonlin[1])) + klin[2])*x[0]
    dx[1] = klin[1]*(1 + x[3]/pow(knonlin[0],knonlin[1]))*x[0] - (klin[3] + klin[4])*x[1]
    dx[2] = klin[3]*x[1] - (klin[5] + klin[6])*x[2]
    dx[3] = klin[5]*x[2] - klin[7]*x[3]
    return(dx)

def func_POSm4_comp(x,klin,knonlin):
    '''
    example function: complex-valued chain model with positive feedback regulation as 
    from Baum et al., 2016
    
    :Parameters:

        - x (numpy array, list, tuple): numerical, ordered values of length 4; these are the 
          variable values.
        - klin: (numpy array, list, tuple): numerical, ordered values of length 8; these are 
          some of the kinetic parameters of the model
        - knonlin: 

    :Returns:

        A complex-valued numpy array of length 4 that contains the time-derivatives of the four variables.

    :Details:

        This function is the same as func_POSm4(), only that the output is complex. Thus, it can
        be used with complex-step derivatives.
        It is a function delivering an ordinary differential equation system with chain structure
        with 4 species and a positive feedback regulation from the last species on the
        conversion between species 1 and 2. The system can give rise to oscillations; 
        see Baum et al., 2016 for details.

    :Reference:
    
        Baum K, Politi AZ, Kofahl B, Steuer R, Wolf J: Feedback, Mass Conservation and Reaction 
        Kinetics Impact the Robustness of Cellular Oscillationss. Plos Comput Biol. 
        2016;12(12):e1005298.
    '''
    dx = np.zeros(4,dtype='complex')
    dx[0] = klin[0]-(klin[1]*(1 + x[3]/pow(knonlin[0],knonlin[1])) + klin[2])*x[0]
    dx[1] = klin[1]*(1 + x[3]/pow(knonlin[0],knonlin[1]))*x[0] - (klin[3] + klin[4])*x[1]
    dx[2] = klin[3]*x[1] - (klin[5] + klin[6])*x[2]
    dx[3] = klin[5]*x[2] - klin[7]*x[3]
    return(dx)



