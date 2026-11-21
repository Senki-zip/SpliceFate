"""Start SpliceFate from YAML; numerical libraries load after thread setup."""
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from splicefate.runner import main

if __name__=='__main__':
    raise SystemExit(main())
