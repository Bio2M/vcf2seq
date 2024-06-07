# __init__.py

import sys, os

app_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, app_path)

from vcf2seq.vcf2seq import compute
