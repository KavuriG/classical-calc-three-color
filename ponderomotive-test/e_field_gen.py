#!usr/bin/python

import numpy as np


ELEC_MASS = 9.10938356E-31
ELEC_CHARGE = -1.60217662E-19
FUND_FREQ = 3.7474057E14
SP_LIGHT = 3E8



class HorizLinPolLight():
    '''creates a linearly polarized light wave with specified params.
    Objects:
        The "optic_dict" object holds the 2 available methods, which makes
        it easy to call these methods programatically. E.g. to call half_wave_plate
        on the beam 'beam_1', you can do: 
                
            beam_1.optic_dict['h'](angle)
            
        instead of
        
             beam_1.half_wave_plate(angle)
    Methods:
        half_wave_plate and quarter_wave_plate apply the relevant phase shifts to the components
        of the beam, and take the angle wrt. horizontal as argument (deg.)

        CAUTION: Setting angle to 0 does NOT make the plate disappear!
                 You'll still see the relevant shifts.
        '''
        
    def half_wave_plate(self, angle=0):
        a_11 = np.cos(np.pi*angle/90)
        a_12= np.sin(np.pi*angle/90)
        a_21 = a_12
        a_22 = -1*a_11
        half_jones_mat = np.exp(-1j*np.pi/2)*np.array([[a_11, a_12],[a_21, a_22]])
        self.jones_vec = np.dot(half_jones_mat, self.jones_vec)
    
    def quarter_wave_plate(self, angle=0):
        sine  = np.sin(np.pi*angle/180)
        cosine = np.cos(np.pi*angle/180)
        a_11 = cosine**2 + (sine**2)*1j
        a_12 = (1 - 1j)*sine*cosine
        a_21 = a_12
        a_22 = sine**2 + (cosine**2)*1j
        quarter_jones_mat = np.exp(-1j*np.pi/4)*np.array([[a_11, a_12],[a_21, a_22]])
        self.jones_vec = np.dot(quarter_jones_mat, self.jones_vec)

    def __init__(self, freq, ampl=1, phase=0., pulse_width=0):
        '''phase in sec.'''
        self.freq = freq
        self.phase = phase
        self.pulse_width = pulse_width
        self.optic_dict = {'h': self.half_wave_plate,
                           'q': self.quarter_wave_plate}
        self.jones_vec = np.array([[ampl*(np.exp(1j*2*np.pi*phase*freq))],[0. +0.j]])



def e_field_gen(num_beams, pulsed=True, *args, **kwargs):
    '''Generates e-field according to specification.

    Args have to be, in order:
        num_beams (say n), pulsed(true), n beam amplitudes, n delays (sec),
        n frequencies, angles[[]], n pulse FWHMs, kwargs of form "b1 = 'hhqh', where
        "1" is beam number and h, q are half and quarter wave plates.'''

    beams = []                   
    for i in range(num_beams):
        beams.append(HorizLinPolLight(args[2*num_beams + i], args[i],
                                      args[num_beams + i],
                                      args[3*num_beams + 1 + i] if pulsed else 0))


    for key, value in kwargs.items():
        '''apply wave plates to beams,
           second argument in kwargs indicates
           beam number. (hacky)'''
        for i in range(len(value)):
            angle = args[3*num_beams][int(key[1]) - 1][i]
            (beams[int(key[1]) - 1].optic_dict[value[i]])(angle)

    if not pulsed:
        '''Add up all the y (or z) components of the beams given, apply envolope
        if pulsed, and return a 2 member list of the functions that can generate
        the field to the required degree of accuracy.'''
        ret = []
        def y_cw(t):
            y_comp = 0
            for beam in beams:
                y_comp += (np.exp(1j*(2*np.pi*(t)*beam.freq - np.pi/4))*beam.jones_vec[0][0])
            #print(beams[0].jones_vec[0][0])
            return np.real(y_comp)
        def z_cw(t):
            z_comp = 0
            for beam in beams:
                z_comp += (np.exp(1j*2*np.pi*(t)*beam.freq)*beam.jones_vec[1][0])
            return  np.real(z_comp)
        ret.append(y_cw)
        ret.append(z_cw)
        return ret

    else:
        ret = []
        def y_pl(t):
            y_comp = 0
            for beam in beams:
                y_comp += np.exp(1j*2*np.pi*(t)*beam.freq)*beam.jones_vec[0][0]*\
                np.exp(-1*2*np.log(2)*((t-beam.phase)/beam.pulse_width)**2)
            return np.real(y_comp)
        def z_pl(t):
#            import pdb; pdb.set_trace()
            z_comp = 0
            for beam in beams:
                z_comp += np.exp(1j*2*np.pi*(t)*beam.freq)*beam.jones_vec[1][0]*\
                np.exp(-1*2*np.log(2)*((t-beam.phase)/beam.pulse_width)**2)
            return np.real(z_comp)
        ret.append(y_pl)
        ret.append(z_pl)
        return ret