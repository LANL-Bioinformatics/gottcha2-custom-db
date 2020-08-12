#!/usr/bin/env python

import sys, os, time
import taxonomy as t

if __name__ == '__main__':
	t.loadTaxonomy( sys.argv[1] if len(sys.argv)>1 else None )

	for line in sys.stdin:
		if line:
			tmp = line.strip().split('\t')
			if not tmp[1]:
				tmp[1] = t.taxid2name(tmp[2])
			print( "\t".join(tmp) )
