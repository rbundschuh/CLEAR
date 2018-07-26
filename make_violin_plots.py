import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

#Import Required Libraries
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

from numpy import sqrt, pi, exp
from lmfit import Model

from pylab import *
from glob import glob

##########################################
###		  From 'fitter.py'			   ###
##########################################

window_size = 250

def dgaussian( x, amp, sig=0.1, cen=0 ):
	out = exp( ( -1.0 * ( x - cen )**2.0 )/( 2.0 * sig**2.0 ) )
	out += exp( ( -1.0 * ( x + cen )**2.0 )/( 2.0 * sig**2.0 ) )
	out *= float( amp ) / float( 2 * sig * sqrt( 2 * pi ) )
	return out

def load_file( name ):
	out = []

	for line in open( name, 'r' ):
		if line.startswith( '#' ) or 'None' in line:
			continue

		line = line.strip().split( '\t' )
		line[ 0 ] = float( line[ 0 ] )
		line[ 1 ] = float( line[ 1 ] )
		line[ 2 ] = int( line[ 2 ] )
		line[ 4 ] = line[ 4 ].startswith( 'T' )
	
		# index, name, count, mu
		#out.append( ( line[ 1 ], line[ 2 ], line[ 3 ] ) )
		#33858.2075218	0.0817306656152	1489	CD74	False
		out.append( ( line[ 3 ], line[ 0 ], line[ 1 ] ) )
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
		model.set_param_hint( 'cen', min=0.0, max=1.0 )
		res = model.fit( bins, x=x, amp=window_size )

		cutoff = 0.5
#		cutoff = 2.5
		if res.params[ 'cen' ].value > cutoff:
			return start

#########################################

font_size = 7

matplotlib.rcParams.update({'font.size': font_size, 'font.family': 'serif'})
matplotlib.rcParams.update({'font.size': font_size, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})

dat_files = glob( '*.dat' )

# Load Profile Tables

profile_tables = []
for f in dat_files:
	to_add = []
	for line in open( str( f ), 'r' ).read().split( '\n' )[1:-1]:
		if 'None' in line: continue
		to_add.append( line.split() )
	profile_tables.append( to_add )

#############################################

to_sort = []

i = 0
for profile in profile_tables:
	cov = 0
	for vector in profile:
		cov += int( vector[ 2 ] )

	to_sort.append( ( cov, profile, dat_files[ i ] ) )
	i += 1

profile_tables = []

i = 0
for a in sorted( to_sort )[::-1]:
	profile_tables.append( a[ 1 ] )
	dat_files[ i ] = a[ 2 ]
	i += 1
del i

#############################################

# Build data vectors
profile_vectors = map(
	lambda x: map(
		lambda y: float( y[1] ), # Select index and ratio
		x
	),
	profile_tables
)

for i in range( len( profile_vectors ) ):
	for j in range( len( profile_vectors[i] ) ):
		profile_vectors[ i ][ j ] = ( j, profile_vectors[ i ][ j ] )

profile_vectors = map( lambda x: map( lambda y: y[1], x ), profile_vectors )

pdf = PdfPages( 'violins.pdf' )

PAGE_WIDTH = 3
PAGE_HEIGHT = 2
page_count = ( len( profile_vectors ) / ( PAGE_WIDTH * PAGE_HEIGHT ) ) + 1

for page_num in range( page_count ):
		fig, axes = plt.subplots(nrows=PAGE_HEIGHT, ncols=PAGE_WIDTH, figsize=(14,7))

		for i in xrange( PAGE_WIDTH * PAGE_HEIGHT ):
				dat_index = page_num * PAGE_WIDTH * PAGE_HEIGHT
				dat_index += i

				x_cord = i % PAGE_WIDTH
				y_cord = i / PAGE_WIDTH

				ax = axes[ y_cord, x_cord ]

				if dat_index >= len( profile_vectors ):
						ax.axis('off')
						ax.get_xaxis().set_visible(False)
						ax.get_yaxis().set_visible(False)
						continue

				dat_file = dat_files[ dat_index ]

				file_reload = load_file( dat_file )
				cut_point = find_state_change( file_reload )
				del file_reload

				dat_file = dat_file.replace( 'dat_files\\', '' )
				dat_file = dat_file.replace( '.dat', '' )
				ax.set_title( dat_file )

				pos, data = [], []
				x_max = 15
				for j in range( 0, x_max ):
						inc = 1000
						p = inc * j
						pos.append( p + (inc / 2) )
						data.append( profile_vectors[ dat_index ][p:p+inc] )

				# Draw the violin plot, note that the bw_method parameter defines the 'definition' of the KDE
				ax.violinplot(data, pos, widths=850, showmeans=True, showextrema=True, showmedians=False, bw_method=0.2)
				del pos, data
				
				if cut_point < x_max * 1000:
					ax.axvline(x=cut_point, color='r', linestyle='-', linewidth=2)

				RIG_FACTOR = 0.1
				ax.set_ylim( -1 - RIG_FACTOR, 1 + RIG_FACTOR )
				ax.set_ylim( -1 - RIG_FACTOR, 1 + RIG_FACTOR )

		fig.tight_layout()
		pdf.savefig( fig )

pdf.close()
