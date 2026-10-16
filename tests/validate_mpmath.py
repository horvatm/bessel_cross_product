"""Independent optional validation: python tests/validate_mpmath.py [build_dir]."""
# SPDX-License-Identifier: GPL-3.0-or-later
import csv
import io
import subprocess
import sys
from pathlib import Path
import mpmath as mp
root = Path(__file__).resolve().parents[1]
build = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else root/'build'
output = subprocess.check_output([str(build/'cross_values')], text=True)
worst = 0
for row in csv.reader(io.StringIO(output)):
    imaginary, v, k, r, actual = map(mp.mpf, row)
    # Imaginary-order products cancel exponentially large imaginary parts.
    mp.mp.dps = max(90, int(2*mp.pi*v/mp.log(10))+60) if imaginary else 100
    order = 1j*v if imaginary else v
    expected = mp.re(mp.besselj(order,k*r)*mp.bessely(order,k)
                     - mp.bessely(order,k*r)*mp.besselj(order,k))
    error = abs(actual-expected)/max(abs(expected),mp.mpf('1e-12'))
    worst = max(worst,float(error))
    print(f"{'imag' if imaginary else 'real':4} v={float(v):7g} k={float(k):7g} "
          f"r={float(r):g} relative_error={float(error):.3g}")
    assert error < mp.mpf('1e-6'), (row, mp.nstr(expected,20))
print(f'{len(output.splitlines())} reference values passed; worst relative error {worst:.3g}')
# Validate example roots independently, including the sign convention.
subprocess.run([str(build/'spectrum'), '20', '.5', '6', str(build/'reference')], check=True)
mp.mp.dps = 90
with (build/'reference.csv').open() as f:
    modes = list(csv.DictReader(f))
for m in modes:
    v=mp.mpf(m['order'])
    imaginary=int(m['imaginary'])
    def boundary(x):
        order=1j*x if imaginary else x
        return mp.re(mp.besselj(order,10)*mp.bessely(order,20)
                     - mp.bessely(order,10)*mp.besselj(order,20))
    reference=mp.findroot(boundary,(v-mp.mpf('.001'),v+mp.mpf('.001')))
    assert abs(v-reference)<mp.mpf('1e-8')
assert sum(int(m['imaginary'])==0 for m in modes)==3
print('All 9 demo orders agree with independent high-precision roots within 1e-8')
