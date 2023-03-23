import os

from dotenv import load_dotenv

"""fealden requies a .env file to specify external file locations and backend to use.

Contents of an example .env file:

FEALDEN_BACKEND=mfold
HYBRID_SS_MIN=/home/username/unafold-new/bin/hybrid-ss-min
SIR_GRAPH=/home/username/mfold/bin/sir_graph
RNASTRUCTURE=/home/username/RNAstructure

FEALDEN_BACKEND can be 'mfold' or 'rnastructure'
"""

load_dotenv()

backend = os.getenv("FEALDEN_BACKEND")

if backend == "mfold":
    from ._unafold import RNAfolder  # noqa
elif backend == "rnastructure":
    from ._rnastructure import RNAfolder  # type: ignore # noqa
else:
    raise ImportError("No backend found, aborting")
