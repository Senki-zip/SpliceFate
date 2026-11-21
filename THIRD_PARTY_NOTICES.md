# Third-party software

The repository's root MIT license applies to the SpliceFate project code. It
does not relicense third-party material under `tools/`.

## CytoTRACE 0.3.1

`tools/CytoTRACE_official/` contains official CytoTRACE source and documentation
from upstream commit `f4ff6c0ab93cc704eef44e17adcd5c6baf09703a`. Source provenance
and checksums are recorded in `COPY_PROVENANCE.json`. SpliceFate invokes the
official R function through a separate R process.

Upstream `DESCRIPTION` declares GPL-3, while upstream `LICENSE.txt` also contains
Stanford non-commercial-use terms. Both files are retained verbatim. Consult
those upstream terms for use or redistribution of CytoTRACE; the SpliceFate MIT
license does not replace them.

Third-party example datasets and upstream test scripts are excluded from this
repository. Other reference checkouts and downloaded executables available in
the development workspace are not included.
