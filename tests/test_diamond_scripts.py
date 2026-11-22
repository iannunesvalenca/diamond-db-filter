"""Unit tests for diamond-db-filter taxonomy utilities and report summarization."""
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import audit_taxids  # noqa: E402
import check_unknown  # noqa: E402
import estimate_distribution  # noqa: E402
import summarize_diamond  # noqa: E402


@pytest.fixture()
def mini_taxdmp(tmp_path):
    """Create a minimal NCBI-format taxonomy dump covering Viruses, Bacteria, Human, and Plant."""
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text(
        "\n".join(
            [
                "1\t|\t1\t|\tno rank\t|\t",
                "131567\t|\t1\t|\tcellular root\t|\t",
                "2\t|\t131567\t|\tdomain\t|\t",
                "561\t|\t2\t|\tgenus\t|\t",
                "562\t|\t561\t|\tspecies\t|\t",
                "2759\t|\t131567\t|\tdomain\t|\t",
                "9605\t|\t2759\t|\tgenus\t|\t",
                "9606\t|\t9605\t|\tspecies\t|\t",
                "33090\t|\t2759\t|\tkingdom\t|\t",
                "3702\t|\t33090\t|\tspecies\t|\t",
                "10239\t|\t1\t|\tsuperkingdom\t|\t",
                "2697049\t|\t10239\t|\tspecies\t|\t",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    names = tmp_path / "names.dmp"
    names.write_text(
        "\n".join(
            [
                "1\t|\troot\t|\t\t|\tscientific name\t|",
                "2\t|\tBacteria\t|\t\t|\tscientific name\t|",
                "561\t|\tEscherichia\t|\t\t|\tscientific name\t|",
                "562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|",
                "2759\t|\tEukaryota\t|\t\t|\tscientific name\t|",
                "9605\t|\tHomo\t|\t\t|\tscientific name\t|",
                "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|",
                "3702\t|\tArabidopsis thaliana\t|\t\t|\tscientific name\t|",
                "10239\t|\tViruses\t|\t\t|\tscientific name\t|",
                "2697049\t|\tSevere acute respiratory syndrome coronavirus 2\t|\t\t|\tscientific name\t|",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    merged = tmp_path / "merged.dmp"
    merged.write_text("999999\t|\t2697049\t|\n", encoding="utf-8")

    return nodes, names, merged


def test_audit_and_unknown_taxids(mini_taxdmp, tmp_path):
    nodes, _, _ = mini_taxdmp
    keep = tmp_path / "keep_taxids.txt"
    keep.write_text("562\n9606\n2697049\n3702\n", encoding="utf-8")

    counts = audit_taxids.audit_taxids(nodes, keep)
    assert counts == {
        "Bacteria": 1,
        "Human": 1,
        "Viruses": 1,
        "Unknown/Other": 1,
    }

    unknown = check_unknown.find_unknown_taxids(nodes, keep)
    assert unknown == [3702]


def test_estimate_distribution_handles_root_self_loop_and_excluded(mini_taxdmp, tmp_path):
    """Ensure TaxIDs descending from root (1 -> 1) terminate cleanly without infinite loops."""
    nodes, _, _ = mini_taxdmp
    keep = tmp_path / "keep_taxids.txt"
    keep.write_text("562\n2697049\n3702\n", encoding="utf-8")

    prot_map = tmp_path / "prot.map"
    prot_map.write_text(
        "\n".join(
            [
                "accession.version\ttaxid",
                "WP_001.1\t562",
                "YP_001.1\t2697049",
                "NP_001.1\t3702",
                "NP_002.1\t9606",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    counts, total_sampled, kept_sampled = estimate_distribution.estimate_distribution(
        nodes, keep, prot_map, sample_rate=1
    )
    assert total_sampled == 4
    assert kept_sampled == 3
    assert counts == {
        "Bacteria": 1,
        "Viruses": 1,
        "Unknown/Other": 1,
        "Excluded": 1,
    }


def test_summarize_diamond_handles_merged_and_blank_taxids(mini_taxdmp, tmp_path):
    nodes, names, merged = mini_taxdmp
    nodes_map = summarize_diamond.load_taxonomy_nodes(nodes)
    names_map = summarize_diamond.load_taxonomy_names(names)
    merged_map = summarize_diamond.load_merged_map(merged)

    tsv = tmp_path / "sample_report.tsv"
    tsv.write_text(
        "\n".join(
            [
                "read1\tYP_001\t98.5\t150\t2\t0\t1\t450\t1\t150\t1e-50\t250.0\t999999;10239\tSARS-CoV-2\tViruses\tViruses",
                "read2\tWP_001\t95.0\t120\t5\t0\t1\t360\t1\t120\t1e-30\t180.0\t562\tE. coli\tBacteria\tBacteria",
                "read3\tXX_001\t80.0\t90\t10\t0\t1\t270\t1\t90\t1e-10\t75.0\t\tUnknown\t\t",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = summarize_diamond.process_diamond_report(
        str(tsv), nodes_map, names_map, merged_map
    )
    assert result is not None
    summary_csv, ranked_csv = result

    summary = pd.read_csv(summary_csv)
    ranked = pd.read_csv(ranked_csv)

    assert len(ranked) == 3
    assert set(summary["Kingdom"]) == {"Viruses", "Bacteria", "Unclassified"}
    # Ensure small e-values are preserved rather than rounded to 0.00
    assert ranked.loc[ranked["qseqid"] == "read1", "evalue"].iloc[0] == pytest.approx(1e-50)
