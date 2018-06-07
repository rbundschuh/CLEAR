import sys
from numpy import sqrt, pi, exp
from lmfit import Model
import matplotlib.pyplot as plt

def dgaussian( x, amp, sig=1.0, cen=2.5 ):
	out = exp( ( -1.0 * ( x - cen )**2.0 )/( 2.0 * sig**2.0 ) )
	out += exp( ( -1.0 * ( x + cen - 10.0 )**2.0 )/( 2.0 * sig**2.0 ) )
	out *= float( amp ) / float( 2 * sig * sqrt( 2 * pi ) )
	return out

def load_file( name ):
	out = []

	for line in open( name, 'r' ):
		if line.startswith( '#' ):
			continue

		line = line.strip().split( '\t' )
		line[ 0 ] = int( line[ 0 ] )
		line[ 2 ] = int( line[ 2 ] )
		line[ 3 ] = float( line[ 3 ] )

		out.append( ( line[ 1 ], line[ 2 ], line[ 3 ] ) )
		del line

	return out

def find_state_change( dataset ):
	x = []
	scale = float( 2 * 10 )
	for i in range( int( scale ) ):
		x.append( i * 10.0 / scale + 0.5 / scale )

	for start in range( len( dataset ) / 10 ):
		start *= 10
	
		a = []
		for i in range( window_size ):
			a.append( dataset[ i + start ][ 2 ] )

		scale = int( scale )
		bins = [ 0 ] * scale
		for i in a:
			bins[ int( i * scale ) ] += 1

		model = Model( dgaussian )
		model.set_param_hint( 'cen', min=0.0, max=5.0 )
		res = model.fit( bins, x=x, amp=window_size )

		cutoff = 2.1
#		cutoff = 2.5
		if res.params[ 'cen' ].value < cutoff:
			return start

data = load_file( sys.argv[ 1 ] )

window_size = 250

for i in range( find_state_change( data ) ):
	print data[i][0]
