from plantvirome.steps.report import run_report


def test_report_builds_from_minimal_outdir(tmp_path):
    d = tmp_path / "s1" / "genbank"
    d.mkdir(parents=True)
    (d / "contigs.qc.tsv").write_text("contig\tlength\tstatus\treasons\nc1\t2700\tsubmit\t\nc2\t900\treject\tshort\n")
    html = run_report(tmp_path).read_text()
    assert "<svg" in html and "s1" in html and "short" in html


def test_report_empty_outdir_does_not_crash(tmp_path):
    assert "Samples" in run_report(tmp_path).read_text()


def test_report_genome_map_with_annotation(tmp_path):
    d = tmp_path / "s1" / "genbank"
    d.mkdir(parents=True)
    (d / "contigs.qc.tsv").write_text("contig\tseqid\tlength\tstatus\treasons\nc1\tc1\t2700\tsubmit\t\n")
    (d / "annotation.tsv").write_text("seqid\torf\tbegin\tend\tstrand\taa\tproduct\tidentity\n"
                                      "c1\t1\t100\t1000\t1\t300\tcoat protein\t99\nc1\t2\t1500\t900\t-1\t200\t\t\n")
    html = run_report(tmp_path).read_text()
    assert "<polygon" in html and "coat protein" in html and "Genome annotation" in html
