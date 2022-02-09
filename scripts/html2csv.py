#!/usr/bin/env python
import sys
import pandas as pd

infile=sys.argv[1]
outfile=sys.argv[2]

with open(infile) as f:
  line=f.readlines()[0]
  t=pd.read_html(line)
  t[0].to_csv(outfile, index=False)
