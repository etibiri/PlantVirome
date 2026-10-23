from plantvirome.steps.report import run_report


def test_report_builds_from_minimal_outdir(tmp_path):
    d = tmp_path / "s1" / "genbank"
    d.mkdir(parents=True)
    (d / "contigs.qc.tsv").write_text("contig\tlength\tstatus\treasons\nc1\t2700\tsubmit\t\nc2\t900\treject\tshort\n")
    html = run_report(tmp_path).read_text()
    assert "<svg" in html and "s1" in html and "short" in html


def test_report_empty_outdir_does_not_crash(tmp_path):
    assert "Samples" in run_report(tmp_path).read_text()
