import sys
from glob import glob

from numpy import sqrt, pi, exp
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import beta as B

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

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
        
        out.append( ( line[ 3 ], line[ 0 ], line[ 1 ] ) )
    
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
            print "Found!", start, "(", res['a'], res['b'], ")"
            return start
    return -1

def sub_file( dataset ):
    return list( ( i[2] for i in dataset ) )

def run_model( bins, x, a, b, amp ):
    model = Model( m )
    
    model.set_param_hint( 'a', min=0, max=15.0, value=a )
    model.set_param_hint( 'b', min=0, max=15.0, value=b )
    model.set_param_hint( 'amp', min=1, max=20, value=amp )
    
    res = model.fit( bins, x=x, nan_policy='omit' )
    
    return res

def run_model( bins, x, a, b, amp ):
    popt, pcov = curve_fit(m, x, bins, bounds=(0, [20,20]))
    return { 'a': popt[0], 'b': popt[1] }#, 'amp': popt[2] }

files = list( sorted( glob( '*.dat' ) ) )
file_data = map( load_file, files )
file_data_sub = map( sub_file, file_data )
limits = map( find_state_change, file_data )

pdf = PdfPages( 'CLEAR_violins.pdf' )

PAGE_WIDTH = 3
PAGE_HEIGHT = 2
page_count = ( len( files ) / ( PAGE_WIDTH * PAGE_HEIGHT ) ) + 1

for page_num in range( page_count ):
	fig, axes = plt.subplots(nrows=PAGE_HEIGHT, ncols=PAGE_WIDTH, figsize=(14,7))

	for i in xrange( PAGE_WIDTH * PAGE_HEIGHT ):
		dat_index = page_num * PAGE_WIDTH * PAGE_HEIGHT +  i

		x_cord = i % PAGE_WIDTH
		y_cord = i / PAGE_WIDTH

		ax = axes[ y_cord, x_cord ]

		if dat_index >= len( files ):
			ax.axis('off')
			ax.get_xaxis().set_visible(False)
			ax.get_yaxis().set_visible(False)
			continue

		ax.set_title( files[dat_index].replace( 'dat_files\\', '' ).replace( '.dat', '' ) )

		pos, data = [], []
		x_max = 15
		for j in range( 0, x_max ):
			inc = 1000
			p = inc * j
			pos.append( p + (inc / 2) )
			a = file_data_sub[ dat_index ][p:p+inc]
			if not len(a): a.append( 0 )
			data.append( a )

		# Draw the violin plot, note that the bw_method parameter defines the 'definition' of the KDE
		ax.violinplot(data, pos, widths=850, showmeans=True, showextrema=True, showmedians=False, bw_method=0.2)
		del pos, data
				
		if limits[ dat_index ] < x_max * 1000:
			ax.axvline( x=limits[ dat_index ], color='r', linestyle='-', linewidth=2 )

		RIG_FACTOR = 0.1
		ax.set_ylim( -1 - RIG_FACTOR, 1 + RIG_FACTOR )
		ax.set_ylim( -1 - RIG_FACTOR, 1 + RIG_FACTOR )

	fig.tight_layout()
	pdf.savefig( fig )

pdf.close()
