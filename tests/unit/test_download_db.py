import gzip
import shutil
import subprocess
import tarfile

import pytest
from simulate import make_tiny_db

from plantvirome.steps.download_db import build_db, discover_parts, gpff_taxids


def test_discover_parts():
    html = ('<a href="viral.1.protein.faa.gz">x</a><a href="viral.1.protein.gpff.gz">y</a>'
            '<a href="viral.wp_protein.2.protein.faa.gz">z</a><a href="viral.1.1.genomic.fna.gz">g</a>')
    assert discover_parts(html) == ["viral.1.protein", "viral.wp_protein.2.protein"]


GPFF = """LOCUS       YP_001             10 aa            linear   VRL 15-OCT-2018
VERSION     YP_001.1
                     /db_xref="taxon:11"
//
LOCUS       YP_002             10 aa            linear   VRL 15-OCT-2018
VERSION     YP_002.1
//
LOCUS       YP_003             10 aa            linear   VRL 15-OCT-2018
VERSION     YP_003.2
                     /db_xref="taxon:12"
//
"""


def test_gpff_taxids(tmp_path):
    f = tmp_path / "x.gpff.gz"
    f.write_bytes(gzip.compress(GPFF.encode()))
    # the record without a taxon xref must not steal the next record's taxid
    assert list(gpff_taxids(f)) == [("YP_001.1", "11"), ("YP_003.2", "12")]


@pytest.mark.skipif(shutil.which("diamond") is None, reason="diamond missing")
def test_build_db_taxonomy_works(tmp_path):
    import random
    from simulate import rand_protein
    rng = random.Random(0)
    p1, p2 = rand_protein(rng), rand_protein(rng)
    src = tmp_path / "src"
    src.mkdir()
    _, tax = make_tiny_db(src, {11: p1, 12: p2})            # only used for its taxdump files
    with tarfile.open(tmp_path / "taxdump.tar.gz", "w:gz") as tf:
        for n in ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"):
            tf.add(tax / n, arcname=n)
    (tmp_path / "a.faa.gz").write_bytes(gzip.compress(f">YP_001.1 d [A]\n{p1}\n>YP_003.2 d [B]\n{p2}\n".encode()))
    (tmp_path / "a.gpff.gz").write_bytes(gzip.compress(GPFF.encode()))
    out = tmp_path / "db"
    dmnd = build_db([tmp_path / "a.faa.gz"], [tmp_path / "a.gpff.gz"], tmp_path / "taxdump.tar.gz", out, 2)
    (tmp_path / "q.faa").write_text(f">q\n{p2}\n")
    res = subprocess.run(["diamond", "blastp", "-d", str(dmnd), "-q", str(tmp_path / "q.faa"),
                          "--outfmt", "6", "qseqid", "staxids"], capture_output=True, text=True, check=True)
    assert res.stdout.split() == ["q", "12"]
