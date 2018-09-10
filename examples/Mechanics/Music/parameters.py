"""Description/definition of the possible setups for
guitar simulation :

- single contact case : one string, one contact. 
- bass guitar : one string, a list of frets.
- fretless guitar : one string, no frets, contacts allowed everywhere on the neck.


Some parameters may be read from input command line (default value depends on the case)
or set by user in notebook :

- fs (sample freq)
- output_freq : frequency of outputs (files)
- restit : restitution coefficient
- final_time : duration of simulation
- matlab_input : where to find case descroption in matlab. Must be "path/prefix"
e.g. for one_contact : 'one_contact/pb1'


"""


### Default parameters of the model (string and contacts) and of the simulation ###

## Case 1 : one contact
one_contact = {'fs': 5000, 'output_freq' : 1, 'restit': 1.,
               'matlab_input': 'one_contact/guitare_obst0', 'nb_modes': 1001, 'length': 1.002,
               'output_name': 'guitar_single_e', 'filt_frets': True,
               'max_coords': (1.5e-3, 0.501)}

## Case 2 : bass guitar, with frets
bass_guitar = {'fs': 15680, 'output_freq' : 1, 'restit': 1.,
               'matlab_input': 'bass_guitar/pb2', 'nb_modes': 862, 'length': .863,
               'output_name': 'bass_e', 'filt_frets': True,
               'max_coords': None}

## Case 3 : fretless bass guitar
fretless_bass_guitar = {'fs': 15680, 'output_freq' : 1, 'restit': 1.,
                        'matlab_input':'fretless_bass_guitar/bsf',
                        'nb_modes': 862, 'length': .863,
                        'output_name': 'fretless_e', 'filt_frets': False,
                        'max_coords': (3.6e-3, 0.64)}

    
    
    
