import argparse

parser = argparse.ArgumentParser(description='Groups transcript lists together.')
parser.add_argument('fnames', metavar='File', type=str, nargs='+',
	help='an integer for the accumulator')
parser.add_argument('--require-samples', metavar='#', type=int, nargs='?',
	help='The number of samples to require a transcript appear in to be included in the output list. Default: number of files given.')

args = parser.parse_args()

all_genes = set()
all_samples = []

print("Processing files...", args.fnames) 
for file in args.fnames:
	with open(file, 'r') as f:
		data = f.read().strip().splitlines()
		for l in data:
			all_genes.add( l )
		all_samples.append( data )

for gene in sorted( all_genes ):
	score = sum( 1 for x in all_samples if gene in x )
	if args.require_samples is not None:
		if score < args.require_samples:
			continue
	else:
		if score is not len(args.fnames):
			continue
	print(gene)