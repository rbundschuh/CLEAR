from collections import defaultdict
import sys
import bisect

if len( sys.argv ) < 3:
	print "Usage: "
	print "\tmake_dat.py [reference txt file] [profiled BED file]"
	sys.exit( 2 )

ref_name = sys.argv[1]
file_name = sys.argv[2]

coverage = defaultdict( lambda: [ [], [], [] ] )
coverage_lengths = defaultdict( lambda: 0 )

for line in open( file_name, 'r' ):
    line = line.strip().split()
    if len( line ) < 4: continue
        
    coverage[ line[0] ][ 0 ].append( int( line[1] ) )
    coverage[ line[0] ][ 1 ].append( int( line[2] ) )
    coverage[ line[0] ][ 2 ].append( int( line[3] ) )

for key in coverage.keys():
    coverage_lengths[ key ] = len( coverage[key][0] )

models = []

for line in open( ref_name, 'r' ):
    line = line.strip().split()
    if len( line ) is 0 or len( line ) is 1: continue
    
    seq_name = line[1]
    chr_name = line[2]
    strand = ('+' in line[3])
    name = line[12]
    exons = []

    starts = map( int, line[9].split(',')[:-1] )
    stops = map( int, line[10].split(',')[:-1] )
    
    if not (len( starts ) == len( stops )):
        print "START/END MISMATCH!!!"
        sys.exit( 0 )
    
    length = 0
    for i in xrange( len( starts ) ):
        exon = ( starts[i], stops[i] )
        length += stops[i] - starts[i]
        exons.append( exon )
        
    toadd = ( seq_name, name, chr_name, strand, exons, length )
    models.append( toadd )

def get_coverage( chromosome, locus ):
    lower = bisect.bisect_left( coverage[ chromosome ][0], locus ) - 1
    #print lower
    try:
        if coverage[ chromosome ][ 0 ][ lower ] <= locus and coverage[ chromosome ][ 1 ][ lower ] > locus:
            return coverage[ chromosome ][ 2 ][ lower ]
        lower += 1
        if coverage[ chromosome ][ 0 ][ lower ] <= locus and coverage[ chromosome ][ 1 ][ lower ] > locus:
            return coverage[ chromosome ][ 2 ][ lower ]
    except:
        return 0
    return 0

def test_gene():
    results = []
    processed = 0
    for model in models:
        i, total, mu = 0, 0, 0
        
        for exon in model[4]:
            for j in xrange( exon[0], exon[1] ):
                c = get_coverage( model[2], j )
                
                total += c
                mu += c * i
                i += 1
        
        if total is not 0:
            mu = 2 * float( mu ) / float( total )
            mu /= float( model[5] )
            mu -= 1
        else:
            mu = None
        
        if mu is not None and not model[3]: #stranded?
            mu *= -1
        
        results.append( ( total, mu, model[5], model[1] ) )
        
        processed += 1
        if processed % 1000 is 0:
            print processed, "processed,", len(models)-processed, "to go! (", 100*float(processed)/len(models), "% )"
    return results

results = sorted( test_gene() )[::-1]

f = open( file_name + '.dat', 'w' )
for a in sorted( results )[::-1]:
    a = map( str, a )
    f.write( "\t".join(a) )
    f.write( "\n" )
f.close()
