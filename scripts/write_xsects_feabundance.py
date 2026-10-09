"""
Write a cross-section file with the iron abundance set to a given value.

The Fe (XKFE*) columns of xsects.dat give the Fe K-shell photoionisation
cross-section (relative line yields are applied on top of it when the file
is loaded); they are scaled by the abundance factor, and the same scaling is
applied to the photoelectric absorption (XPHOT) above the Fe K edge, where
Fe dominates the absorption.

The written file can be used by pointing the XARS_XSECTS_FILE environment
variable to it, e.g.:

  python scripts/write_xsects_feabundance.py --abundance 10 --output xsects_fe10.dat
  XARS_XSECTS_FILE=xsects_fe10.dat python examples/sphere_feabundance.py
"""

import argparse
import os

import numpy

parser = argparse.ArgumentParser(
    description="""Write a cross-section file with a modified iron abundance""")
parser.add_argument('--abundance', type=float, required=True,
                    help='Fe abundance relative to solar (1 = solar)')
parser.add_argument('--input', type=str, default=None,
                    help='input cross-section file (default: xsects.dat shipped with xars)')
parser.add_argument('--output', type=str, default=None,
                    help='output file (default: xsects_fe<abundance>.dat)')
args = parser.parse_args()

input_file = args.input
if input_file is None:
    import xars.xsects
    input_file = os.path.join(os.path.dirname(xars.xsects.__file__), 'xsects.dat')
output_file = args.output
if output_file is None:
    output_file = 'xsects_fe%g.dat' % args.abundance

with open(input_file) as f:
    lines = f.readlines()

# column names from the table header comment
# (e.g. "# E XPHOT XKFEa2 XKFEa1 XKFEb XKC XKO ... XKNIa2 XKNIa1")
colnames = [l for l in lines if l.startswith('# E XPHOT')][-1].lstrip('#').split()
fe_columns = numpy.array(
    [2 + i for i, name in enumerate(colnames[2:]) if name.startswith('XKFE')],
    dtype=int)
assert len(fe_columns) > 0, 'no iron (XKFE*) columns found in %s' % input_file

# the table follows the comments and the line energy/yield/width/asymmetry rows
data_rows = [i for i, l in enumerate(lines) if not l.startswith('#')]
nheader_rows = 4
table_start = data_rows[nheader_rows]

table = numpy.loadtxt(input_file)
z_fe = args.abundance
fe_photo = table[nheader_rows:, fe_columns[0]].copy()
table[nheader_rows:, fe_columns] *= z_fe
table[nheader_rows:, 1] += (z_fe - 1) * fe_photo
assert (table >= 0).all(), 'scaling made the cross-sections negative'

with open(output_file, 'w') as f:
    f.write(''.join(lines[:table_start]))
    numpy.savetxt(f, table[nheader_rows:])
print('wrote %s with Fe abundance %g (from %s)' % (output_file, z_fe, input_file))
