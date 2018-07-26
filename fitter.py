import sys
from glob import glob

from numpy import sqrt, pi, exp
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import beta as B

window_size = 250

def m( x, a, b ):#, amp ):
    #if x < 0 or x > 1: return 0
    x = ( x + 1.0 ) / 2.0
    
    out = ( x * (1 - x)**(1 + b) ) / B(2,2 + b) 
    out += ( (x)**(1 + a) ) / B(2 + a,2) 
    out -= ( (x)**(2 + a) ) / B(2 + a,2)
    amp = (window_size * 0.1) / 4.0 # Normalization Factor
    out *= amp
    
    if type( x ) is np.ndarray or type( x ) is np.array:
        for i in xrange( len( x ) ):
            if x[i] < -1 or x[i] > 1:
                out[i] = 0
    elif (type( x ) is int or type(x) is float) and ( x < -1 or x > 1 ):
        out = 0
    return out 

def load_file( name ):
    out = []

    for line in open( name, 'r' ):
        if line.startswith( '#' ) or 'None' in line: continue

        line = line.strip().split( '\t' )
        line[ 0 ] = float( line[ 0 ] )
        line[ 1 ] = float( line[ 1 ] )
        line[ 2 ] = int( line[ 2 ] )
        line[ 4 ] = line[ 4 ].startswith( 'T' )

        # index, name, count, mu
        #out.append( ( line[ 1 ], line[ 2 ], line[ 3 ] ) )
        #33858.2075218	0.0817306656152	1489	CD74	False
        out.append( ( line[ 3 ], line[ 0 ], line[ 1 ] ) )
        #out.append( ( line[ 1 ] )

    return out

def find_state_change( dataset ):
    for start in range( len( dataset ) / 10 ):
        start *= 10

        a = [ i[2] for i in dataset[ start:start+window_size] ]
        bin_count = 20
        
        bins = np.histogram( a, bins=bin_count, range=[-1,1] )
        x = bins[1][:-1] + (1.0/bin_count)
        bins = bins[0] 
        
        res = run_model( bins, x, 0, 0, 5 )
        if res[ 'a' ] > 2 or res[ 'b' ] > 2: 
            #print "Transition Found!", start, "(", res['a'], res['b'], ")"
            return start
    return -1

def sub_file( dataset ):
    return list( ( i[2] for i in dataset ) )

def run_model( bins, x, a, b, amp ):
    popt, pcov = curve_fit(m, x, bins, bounds=(0, [20,20]))
    return { 'a': popt[0], 'b': popt[1] }#, 'amp': popt[2] }

data = load_file( sys.argv[ 1 ] )
change = find_state_change( data )

#print change 
data = [ i[0] for i in data ]
print '\n'.join( data[:change] )
