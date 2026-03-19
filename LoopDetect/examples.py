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
        - knonlin: (numpy array, list, tuple): numerical, ordered values of length 2; these are 
          some of the non-linear kinetic parameters of the model

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
    fb = 1 + (x[3] / knonlin[0])**knonlin[1]
    dx = np.zeros(4)
    dx[0] = klin[0] - klin[1] * x[0] * fb - klin[2] * x[0]
    dx[1] = klin[1] * x[0] * fb - klin[3] * x[1] - klin[4] * x[1]
    dx[2] = klin[3] * x[1] - klin[5] * x[2] - klin[6] * x[2]
    dx[3] = klin[5] * x[2] - klin[7] * x[3]
    return dx



def func_POSm4_comp(x,klin,knonlin):
    '''
    example function: complex-valued chain model with positive feedback regulation as 
    from Baum et al., 2016
    
    :Parameters:

        - x (numpy array, list, tuple): numerical, ordered values of length 4; these are the 
          variable values.
        - klin: (numpy array, list, tuple): numerical, ordered values of length 8; these are 
          some of the kinetic parameters of the model
        - knonlin: (numpy array, list, tuple): numerical, ordered values of length 2; these are 
          some of the non-linear kinetic parameters of the model


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
    fb = 1 + (x[3] / knonlin[0])**knonlin[1]
    dx = np.zeros(4)
    dx[0] = klin[0] - klin[1] * x[0] * fb - klin[2] * x[0]
    dx[1] = klin[1] * x[0] * fb - klin[3] * x[1] - klin[4] * x[1]
    dx[2] = klin[3] * x[1] - klin[5] * x[2] - klin[6] * x[2]
    dx[3] = klin[5] * x[2] - klin[7] * x[3]
    return dx.astype(complex)

def func_complex_formation(y, params):
    A, B, AB = y
    k1, k2 = params
    dAdt = -k1 * A * B + k2 * AB
    dBdt = -k1 * A * B + k2 * AB
    dABdt = k1 * A * B - k2 * AB
    
    return np.array([dAdt, dBdt, dABdt], dtype=float)


def markevich04(y, params):
    '''example function: complex formation with positive feedback regulation as from Markevich et al., 2004
    
    :Parameters:
        - y (numpy array, list, tuple): numerical, ordered values of length 4; these are the variable values.
        - params: (numpy array, list, tuple): numerical, ordered values of length 7; these are the kinetic parameters of the model
        
    :Returns:
        A numpy array of length 4 that contains the time-derivatives of the four variables.
    :Details:
        This function encodes an ordinary differential equation system with 4 species and a positive feedback regulation from the last species on the conversion between species 1 and 2. The system can give rise to oscillations; see Markevich et al., 2004 for details.
    :Reference:
        Markevich NI, Hoek JB, Kholodenko BN: Signaling switches and bistability arising from multisite phosphorylation in protein kinase cascades. J Cell Biol. 2004;164(3):353-359.
        Varusai TM, Kolch W, Kholodenko BN, Nguyen LK: Protein-protein interactions generate hidden feedback and feed-forward loops to trigger bistable switches, oscillations and biphasic dose-responses. Mol BioSyst. 2015;11(10):2750-2762.
    '''
    A = y[0]
    B = y[1]
    AB = y[2]
    B_star = y[3]
    k1, k2, k3, k4, kn1, kn2, kn3 = params
    dAdt = -k1 * A * B + k2 * AB
    dBdt = -k1 * A * B + k2 * AB + k3 * A * (B_star / (B_star + kn1)) - k4 * kn3 * (B / (B + kn2))
    dABdt = k1 * A * B - k2 * AB
    dB_stardt = -k3 * A * (B_star / (B_star + kn1)) + k4 * kn3 * (B / (B + kn2))
    return np.array([dAdt, dBdt, dABdt, dB_stardt], dtype=float)

def func_huang96_mapk(y, k):
    '''
    example function: complex formation with positive feedback regulation as from Huang and Ferrell, 1996
    :Parameters:
        - y (numpy array, list, tuple): numerical, ordered values of length 22; these are the variable values.
        - k: (numpy array, list, tuple): numerical, ordered values of length 30; these are the kinetic parameters of the model
    :Returns:
        A numpy array of length 22 that contains the time-derivatives of the twenty-two variables.
    :Details:
        This function encodes an ordinary differential equation system with 22 species, to model the mitogen-activated protein kinase cascade; see Huang and Ferrell, 1996 for details.
    :Reference:
        Huang W, Ferrell JE: Ultrasensitivity and bistability in the mitogen-activated protein kinase cascade. Proc Natl Acad Sci U S A. 1996;93(14):7240-7245.
    '''
    S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, \
    S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, \
    S21, S22 = y

    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, \
    k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, \
    k21, k22, k23, k24, k25, k26, k27, k28, k29, k30 = k

    dS1  = -k1*S1*S3 + k2*S13 + k3*S13
    dS2  = -k4*S2*S4 + k5*S14 + k6*S14
    dS3  = -k1*S1*S3 + k2*S13 + k6*S14

    dS4  = (k3*S13 - k4*S2*S4 + k5*S14
            - k7*S5*S4 + k8*S15 + k9*S15
            - k13*S6*S4 + k14*S16 + k15*S16)

    dS5  = -k7*S5*S4 + k8*S15 + k12*S20

    dS6  = (k9*S15 - k10*S6*S12 + k11*S20
            - k13*S6*S4 + k14*S16 + k18*S19)

    dS7  = (k15*S16 - k16*S7*S12 + k17*S19
            - k19*S8*S7 + k20*S17 + k21*S17
            - k25*S9*S7 + k26*S18 + k27*S18)

    dS8  = -k19*S8*S7 + k20*S17 + k24*S22

    dS9  = (k21*S17 - k22*S9*S11 + k23*S22
            - k25*S9*S7 + k26*S18 + k30*S21)

    dS10 = k27*S18 - k28*S10*S11 + k29*S21

    dS11 = (-k22*S9*S11 + k23*S22 + k24*S22
            - k28*S10*S11 + k29*S21 + k30*S21)

    dS12 = (-k10*S6*S12 + k11*S20 + k12*S20
            - k16*S7*S12 + k17*S19 + k18*S19)

    dS13 = k1*S1*S3 - k2*S13 - k3*S13
    dS14 = k4*S2*S4 - k5*S14 - k6*S14
    dS15 = k7*S5*S4 - k8*S15 - k9*S15
    dS16 = k13*S6*S4 - k14*S16 - k15*S16
    dS17 = k19*S8*S7 - k20*S17 - k21*S17
    dS18 = k25*S9*S7 - k26*S18 - k27*S18
    dS19 = k16*S7*S12 - k17*S19 - k18*S19
    dS20 = k10*S6*S12 - k11*S20 - k12*S20
    dS21 = k28*S10*S11 - k29*S21 - k30*S21
    dS22 = k22*S9*S11 - k23*S22 - k24*S22

    return np.array([
        dS1, dS2, dS3, dS4, dS5, dS6, dS7, dS8, dS9, dS10,
        dS11, dS12, dS13, dS14, dS15, dS16, dS17, dS18,
        dS19, dS20, dS21, dS22
    ], dtype=float)

def func_kho00(y, klin, kn, n1):
    '''example function: complex formation with positive feedback regulation as from Khoo et al., 2000
    :Parameters:
        - y (numpy array, list, tuple): numerical, ordered values of length 8; these are the variable values.
        - klin: (numpy array, list, tuple): numerical, ordered values of length 10; these are some of the kinetic parameters of the model
        - kn: (numpy array, list, tuple): numerical, ordered values of length 11; these are some of the saturation constants of the model
    :Returns:
        A numpy array of length 8 that contains the time-derivatives of the eight variables.
    :Details:
        The model for the mitogen-activated phosphorylation kinase (MAPK) cascade proposed in [Reference] captures the first, second and 
        third level of a phosphorylation cascade which occurs for example in epidermal growth factor (EGF) signaling. The kinase on the 
        first level can be reversibly phosphorylated once, the kinases on levels two andthree can be reversibly phosphorylated twice. 
        The phosphorylated kinase of the first level acts positively on the phosphorylations on the second level, the double phosphorylated 
        kinase of the second level acts positively on the phosphorylations on the third level. The model contains an explicit negative feedback regulation: 
        The double phosphorylated kinase of the third level acts negatively on the phosphorylation on the first level. All phosphorylations 
        are governed by Michaelis-Menten kinetics. The model contains eight variables, ten rate coefficients, eleven inhibition or 
        activation constants, and one Hill coefficient.
    :Reference:
        Khoo MS, McCluskey A, McCluskey S, Houslay MD, Baillie GS: Oscillating cyclic AMP levels in cells are regulated by phosphodiesterases and modulated by their phosphorylation status. J Biol Chem. 2000;275(44):33899-33904.
    '''

    S1, S2, S3, S4, S5, S6, S7, S8 = y
    
    # Unpack kinetic parameters
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = klin
    
    # Unpack saturation constants
    kn1, kn2, kn3, kn4, kn5, kn6, kn7, kn8, kn9, kn10, kn11 = kn
    
    # dS1/dt
    term_inhibit = (kn1**n1) / (kn1**n1 + S8**n1)
    dS1 = -k1 * (S1/(kn2 + S1)) * term_inhibit \
          + k2 * (S2/(S2 + kn3))
    
    # dS2/dt
    dS2 =  k1 * (S1/(kn2 + S1)) * term_inhibit \
          - k2 * (S2/(S2 + kn3))
    
    # dS3/dt
    dS3 = k6 * (S4/(kn7 + S4)) \
          - k3 * (S2*S3/(kn4 + S3))
    
    # dS4/dt
    dS4 = k3 * (S2*S3/(kn4 + S3)) \
          - k6 * (S4/(kn7 + S4)) \
          + k5 * (S5/(kn6 + S5)) \
          - k4 * (S2*S4/(kn5 + S4))
    
    # dS5/dt
    dS5 = k4 * (S2*S4/(kn5 + S4)) \
          - k5 * (S5/(kn6 + S5))
    
    # dS6/dt
    dS6 = k10 * (S7/(kn11 + S7)) \
          - k7 * (S5*S6/(kn8 + S6))
    
    # dS7/dt
    dS7 =  k7 * (S5*S6/(kn8 + S6)) \
            - k10 * (S7/(kn11 + S7)) \
            - k8 * (S5*S7/(kn9 + S7)) \
            + k9 * (S8/(kn10 + S8))
    
    # dS8/dt
    dS8 =  k8 * (S5*S7/(kn9 + S7)) \
           - k9 * (S8/(kn10 + S8))
    
    return np.array([dS1, dS2, dS3, dS4, dS5, dS6, dS7, dS8], dtype=float)

def func_kho00_nofb(y, klin, kn, n1):
    '''example function: complex formation without feedback regulation as from Khoo et al., 2000
    :Parameters:
        - y (numpy array, list, tuple): numerical, ordered values of length 8; these are the variable values.
        - klin: (numpy array, list, tuple): numerical, ordered values of length 10; these are some of the kinetic parameters of the model
        - kn: (numpy array, list, tuple): numerical, ordered values of length 11; these are some of the saturation constants of the model
        - n1: (float): Hill coefficient for the inhibition of the phosphorylation on the first level by the double phosphorylated kinase of the third level 
          not used in this function, as there is no feedback regulation included, but required for func_kho00())
    :Returns:
        A numpy array of length 8 that contains the time-derivatives of the eight variables.
    :Details:
        This function encodes a model for the mitogen-activated phosphorylation kinase (MAPK) cascade proposed in [Reference] that captures the first, second and third level of a phosphorylation cascade which occurs for example in epidermal growth factor (EGF) signaling. The kinase on the first level can be reversibly phosphorylated once, the kinases on levels two and three can be reversibly phosphorylated twice. The phosphorylated kinase of the first level acts positively on the phosphorylations on the second level, the double phosphorylated kinase of the second level acts positively on the phosphorylations on the third level. The model contains an explicit negative feedback regulation: The double phosphorylated kinase of the third level acts negatively on the phosphorylation on the first level. In contrast to func_khoo00(), this function does not contain this negative feedback regulation. All phosphorylations are governed by Michaelis-Menten kinetics. The model contains eight variables, ten rate coefficients and eleven inhibition or activation constants.
    :Reference:
        Khoo MS, McCluskey A, McCluskey S, Houslay MD, Baillie GS: Oscillating cyclic AMP levels in cells are regulated by phosphodiesterases and modulated by their phosphorylation status. J Biol Chem. 2000;275(44):33899-33904.
    '''
    
    S1, S2, S3, S4, S5, S6, S7, S8 = y
    
    # Unpack kinetic parameters
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = klin
    
    # Unpack saturation constants
    kn1, kn2, kn3, kn4, kn5, kn6, kn7, kn8, kn9, kn10, kn11 = kn
    
    # dS1/dt
    # term_inhibit = (kn1**n1) / (kn1**n1 + S8**n1)
    dS1 = -k1 * (S1/(kn2 + S1)) \
          + k2 * (S2/(S2 + kn3))
    
    # dS2/dt
    dS2 =  k1 * (S1/(kn2 + S1)) \
          - k2 * (S2/(S2 + kn3))
    
    # dS3/dt
    dS3 = k6 * (S4/(kn7 + S4)) \
          - k3 * (S2*S3/(kn4 + S3))
    
    # dS4/dt
    dS4 = k3 * (S2*S3/(kn4 + S3)) \
          - k6 * (S4/(kn7 + S4)) \
          + k5 * (S5/(kn6 + S5)) \
          - k4 * (S2*S4/(kn5 + S4))
    
    # dS5/dt
    dS5 = k4 * (S2*S4/(kn5 + S4)) \
          - k5 * (S5/(kn6 + S5))
    
    # dS6/dt
    dS6 = k10 * (S7/(kn11 + S7)) \
          - k7 * (S5*S6/(kn8 + S6))
    
    # dS7/dt
    dS7 =  k7 * (S5*S6/(kn8 + S6)) \
            - k10 * (S7/(kn11 + S7)) \
            - k8 * (S5*S7/(kn9 + S7)) \
            + k9 * (S8/(kn10 + S8))
    
    # dS8/dt
    dS8 =  k8 * (S5*S7/(kn9 + S7)) \
           - k9 * (S8/(kn10 + S8))
    
    return np.array([dS1, dS2, dS3, dS4, dS5, dS6, dS7, dS8], dtype=float)

def func_gold90(y, klin, kn, n):
    '''
    example function: calcium oscillations with positive feedback regulation as from Goldbeter, 1990
    :Parameters:
        - y (numpy array, list, tuple): numerical, ordered values of length 2; these are the variable values.
        - klin: (numpy array, list, tuple): numerical, ordered values of length 6; these are some of the kinetic parameters of the model
        - kn: (numpy array, list, tuple): numerical, ordered values of length 3; these are some of the saturation constants of the model
        - n: (numpy array, list, tuple): numerical, ordered values of length 3; these are some of the Hill coefficients of the model
    :Returns:
        A numpy array of length 2 that contains the time-derivatives of the two variables.
    :Details:
        The model for calcium oscillations proposed in [Reference] captures the dynamics of calcium concentration in the cytosol and in the endoplasmic reticulum (ER). 
        Calcium is released from the ER into the cytosol by a calcium channel that is activated by calcium itself, thus giving rise to a positive feedback regulation. 
        Calcium is pumped back into the ER by a calcium pump. The model contains two variables, six kinetic parameters, three saturation constants and three Hill coefficients.
    :Reference: 
        Goldbeter, A. (1990). A model for circadian rhythms based on the dynamics of calcium oscillations in the cytosol. Proceedings of the National Academy of Sciences, 87(19), 7264-7268.
    '''
    S1, S2 = y
    k1, k2, k3, k4, k5, k6 = klin
    kn1, kn2, kn3 = kn
    n1, n2, n3 = n

    v1 = k1
    v2 = k2
    v3 = k3 * (S1**n1) / (kn1**n1 + S1**n1)
    v4 = k4 * (S2**n2)/(kn2**n2 + S2**n2) * (S1**n3)/(kn3**n3 + S1**n3)
    v5 = k5 * S2
    v6 = k6 * S1

    dS1 = v1 + v2 - v3 + v4 + v5 - v6
    dS2 = v3 - v4 - v5

    return np.array([dS1, dS2])

def func_elo00(y, klin, knonlin):
    '''
    example function: gene regulatory network with positive feedback regulation as from Elowitz and Leibler, 2000 (Repressilator)
    :Parameters:
        - y (numpy array, list, tuple): numerical, ordered values of length 6; these are the variable values.
        - klin: (numpy array, list, tuple): numerical, ordered values of length 15; these are some of the kinetic parameters of the model
        - knonlin: (numpy array, list, tuple): numerical, ordered values of length 6; these are some of the non-linear kinetic parameters of the model
    :Returns:
        A numpy array of length 6 that contains the time-derivatives of the six variables.
    :Details:
        The model for a synthetic gene regulatory network proposed in [Reference] captures the dynamics of three genes that repress each other in a cyclic manner. 
        Each gene is transcribed into mRNA and then translated into protein. The transcription of each gene is inhibited by the protein product of another gene, thus giving rise to a positive feedback regulation. 
        The model contains six variables (three for the proteins and three for the mRNAs) and fifteen kinetic parameters.
    :Reference:
        Elowitz MB, Leibler S: A synthetic oscillatory network of transcriptional regulators. Nature. 2000;403(6767):335-338.
    '''
    LacI, TetR, lambda_cl, LacI_mRNA, TetR_mRNA, lambda_cl_mRNA = y
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15 = klin
    kn1, kn2, kn3, kn4, kn5, kn6 = knonlin

    dLacI_dt = -k1 * LacI + k2 * LacI_mRNA
    dTetR_dt = -k3 * TetR + k4 * TetR_mRNA
    dlambda_cl_dt = -k5 * lambda_cl + k6 * lambda_cl_mRNA

    dLacI_mRNA_dt = k7 + k8 * (kn1**kn4) / (kn1**kn4 + lambda_cl**kn4) - k9 * LacI_mRNA
    dTetR_mRNA_dt = k10 + k11 * (kn2**kn5) / (kn2**kn5 + LacI**kn5) - k12 * TetR_mRNA
    dlambda_cl_mRNA_dt = k13 + k14 * (kn3**kn6) / (kn3**kn6 + TetR**kn6) - k15 * lambda_cl_mRNA

    return np.array([dLacI_dt, dTetR_dt, dlambda_cl_dt, dLacI_mRNA_dt, dTetR_mRNA_dt, dlambda_cl_mRNA_dt])