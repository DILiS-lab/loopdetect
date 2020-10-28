import pkg_resources
import numpy as np
import pandas as pd

def func_li08(x,t):
    '''example function: bacterial cell cycle [modelwtin(t,y), Li et al. 2008]
    
    :Parameters:

        - x (numpy array, list, tuple): numerical of length 18, ordered input of the variable values
        - t (float): time variable. Note that the function is independent from time.

    :Returns:

        A numpy array of length 18 containing the temporal derivative of the 18 variables.

    :Details:
        
        This is a function encoding an ordinary differential equation model
        that delivers the dynamics of the 
        Caulobacter cell cycle. Note that to obtain the solution as published, also events have to be 
        considered, i.e. certain conditions lead to a change in certain variable 
        values; see Li et al., 2008 for details.
    
    :Reference:

        Li S, Brazhnik P, Sobral B, Tyson JJ. A Quantitative Study of the 
        Division Cycle of Caulobacter crescentus Stalked Cells. Plos Comput Biol. 
        2008;4(1):e9.
    '''
    #
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
    #Parameters values for the equations
    #
    ksCtrAP1  = 0.0083 
    JiCtrACtrA  = 0.4 	
    niCtrACtrA  = 2 
    ksCtrAP2  = 0.073 	
    JaCtrACtrA  = 0.45 	
    naCtrACtrA  = 2 
    kdCtrA1  = 0.002 
    kdCtrA2  = 0.15 	
    ndCtrA2  = 2 
    JdCtrADivKP  = 0.55 
    ksGcrA  = 0.045 
    JiGcrACtrA  = 0.2 	
    niGcrACtrA  = 2 
    kdGcrA  = 0.022 
    ksFts  = 0.063 	
    kdFts  = 0.035 
    kzringopen  = 0.8 	
    Jaopen  = 0.01 
    kzringclosed1  = 0.0001 	
    Jaclosed1  = 0.1 
    kzringclosed2  = 0.6 	
    nzringclosed2  = 4 
    JZringFts  = 0.78 
    ksDivK  = 0.0054 	
    ktransDivKP  = 0.0295 	
    ktransDivK  = 0.5 	
    kdDivK  = 0.002 
    ksI  = 0.08 	
    kdI  = 0.04 
    ksCcrM  = 0.072 	
    kdCcrM  = 0.07 
    kaDnaA  = 0.0165 	
    JiDnaAGcrA  = 0.5 	
    niDnaAGcrA  = 2 
    kdDnaA  = 0.007 
    kaIni  = 0.01 	
    JaIni  = 1 		
    naIni  = 4 
    thetaCtrA  = 0.2 	
    nthetaCtrA  = 4 
    thetaDnaA  = 0.6 	
    nthetaDnaA  = 4 
    thetaGcrA  = 0.45 	
    nthetaGcrA  = 4 
    thetaCori  = 0.0002 	
    nthetaCori  = 1 
    kmcori  = 0.4 	
    Jmcori  = 0.95 	
    nmcori  = 4 
    kmccrM  = 0.4 	
    JmccrM  = 0.95 	
    nmccrM  = 4 
    kmctrA  = 0.4 	
    JmctrA  = 0.95 	
    nmctrA  = 4 
    kmfts  = 0.4 	
    Jmfts  = 0.95 	
    nmfts  = 4 
    kelong  = 0.95/160 
    nelong  = 4 
    #
    #$end of parameters
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
    ##Differential equations of the model
    ##variable map
    #y[0]  = [CtrA]
    #y[1]  = [GcrA]
    #y[2]  = [DnaA]
    #y[3]  = [Fts]
    #y[4]  = [Zring]
    #y[5]  = [DivK]
    #y[6]  = [DivK~P]
    #y[7]  = [Total DivK]
    #y[8]  = [I] (Intermediate)]
    #y[9]  = [CcrM]
    #y[10]  = [hCori] (hemimethylated Cori)
    #y[11]  = [hccrM] (hemimethylated ccrM)
    #y[12]  = [hctrA] (hemimethylated ctrA)
    #y[13]  = [hfts](hemimethylated fts)
    #y[14]  = [Ini] (Initiation)
    #y[15]  = [Elong](Elongation)
    #y[16]  = [DNA] (Total DNA)
    #y[17]  = Count (# of Chromosome) 
    #
    dydt  = np.zeros(18)
    y=x
    dydt[0]  = (ksCtrAP1*pow(JiCtrACtrA,niCtrACtrA)/(pow(JiCtrACtrA,niCtrACtrA)+pow(y[0],niCtrACtrA))*y[1]+ksCtrAP2*pow(y[0],naCtrACtrA)/(pow(JaCtrACtrA,naCtrACtrA)+pow(y[0],naCtrACtrA)))*y[12]-(kdCtrA1+kdCtrA2*pow(y[6],ndCtrA2)/(pow(JdCtrADivKP,ndCtrA2)+pow(y[6],ndCtrA2)))*y[0] 
    dydt[1]  = (ksGcrA*pow(JiGcrACtrA,niGcrACtrA)/(pow(JiGcrACtrA,niGcrACtrA)+pow(y[0],niGcrACtrA))*y[2]-kdGcrA*y[1]) 
    dydt[2]  = kaDnaA*pow(JiDnaAGcrA,niDnaAGcrA)/(pow(JiDnaAGcrA,niDnaAGcrA)+pow(y[1],niDnaAGcrA))*y[0]*(2-y[10])-kdDnaA*y[2] 
    dydt[3]  = ksFts*y[0]*y[13]-kdFts*y[3] 
    dydt[4]  = (kzringopen*(1-y[4])/(0.01+(1-y[4]))-(kzringclosed1+kzringclosed2*pow((y[3]/JZringFts),nzringclosed2))*y[4]/(0.05+y[4])) 
    dydt[5]  = (ksDivK*y[0]+ktransDivKP*y[6]-ktransDivK*(1-y[4])*y[5]-kdDivK*y[5]) 
    dydt[6]  = (-ktransDivKP*y[6]+ktransDivK*(1-y[4])*y[5]-kdDivK*y[6]) 
    dydt[7]  = (ksDivK*y[0]-kdDivK*y[7]) 
    dydt[8]  = ksI*y[11]*y[0]-kdI*y[8] 
    dydt[9]  = ksCcrM*y[8]-kdCcrM*y[9] 
    dydt[10]  = -kmcori*pow(y[9],nmcori)/(pow(Jmcori,nmcori)+pow(y[9],nmcori))*y[10] 
    dydt[11]  = -kmccrM*pow(y[9],nmccrM)/(pow(JmccrM,nmccrM)+pow(y[9],nmccrM))*y[11] 
    dydt[12]  = -kmctrA*pow(y[9],nmctrA)/(pow(JmctrA,nmctrA)+pow(y[9],nmctrA))*y[12] 
    dydt[13]  = -kmfts*pow(y[9],nmfts)/(pow(Jmfts,nmfts)+pow(y[9],nmfts))*y[13] 
    dydt[14]  = kaIni*pow((y[2]/thetaDnaA),nthetaDnaA)*pow((y[1]/thetaGcrA),4)/(pow(JaIni,naIni)+pow((y[0]/thetaCtrA),nthetaCtrA)+pow((y[2]/thetaDnaA),nthetaDnaA)+pow((y[1]/thetaGcrA),nthetaGcrA)+pow((y[10]/thetaCori),nthetaCori)) 
    dydt[15]  = kelong*pow(y[15],nelong)/(pow(y[15],nelong)+pow(0.05,nelong))*y[17] 
    dydt[16]  = kelong*pow(y[15],nelong)/(pow(y[15],nelong)+pow(0.05,nelong))*y[17] 
    dydt[17]  = 0  #count (of chromosome - is only altered at an event
    #
    return(dydt)
    ##end of equations
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



def load_li08_sol():
    """Loads solutions for an example function of the bacterial cell cycle, Li et al. 2008, 
    as dataframe.
    
    :Returns: 
        A pandas dataframe that provides the solution for 18 species over time (roughly 3 cycles, 
        up to 360 minutes, 634 time points in between). It contains the following columns: 
            
            - time, y1=x0, y2=x1,..., y18=x17

        All are numbers. Time is given in minutes. No input argument required.

    :Reference:

        Li S, Brazhnik P, Sobral B, Tyson JJ. A Quantitative Study of the 
        Division Cycle of Caulobacter crescentus Stalked Cells. Plos Comput Biol. 
        2008;4(1):e9.

    :Example:
        
    Load the content into a variable::
        
        #import the package
        import loopdetect.examples 
        #load the data
        sol_vec = loopdetect.examples.load_li08_sol()

    """
    stream = pkg_resources.resource_stream('loopdetect', 'data/li08_solution.tsv')
    return pd.read_csv(stream, sep='\t', encoding='latin-1')


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




