import pandas as pd

from rdna_miner.steps.compile_taxonomy import abundance_domains, choose_taxonomy


def _pr2(query, domain, confidence, taxon_long, genus="Pr2Genus"):
    return {
        "query": query,
        "domain": domain,
        "domain_confidence": confidence,
        "taxon_long": taxon_long,
        "confidence_long": f"100;{confidence}",
        "class": "Pr2Class",
        "genus": genus,
        "supergroup": "Pr2Supergroup",
    }


def _silva(query, domain, confidence, taxon_long, genus="SilvaGenus"):
    return {
        "query": query,
        "domain": domain,
        "domain_confidence": confidence,
        "taxon_long": taxon_long,
        "confidence_long": f"100;{confidence}",
        "class": "SilvaClass",
        "genus": genus,
        "phylum": "SilvaPhylum",
    }


def _frame(rows):
    return pd.DataFrame(rows)


def _chosen(result, query):
    return result.set_index("query").loc[query]


def test_agreed_bacteria_keeps_silva_even_when_pr2_is_more_confident():
    result = choose_taxonomy(
        _frame([_pr2("bac1", "Bacteria", 90, "Root;Bacteria;Pr2;Genus")]),
        _frame([_silva("bac1", "Bacteria", 40, "Root;Bacteria;Pseudomonadota;Genus")]),
    )

    chosen = _chosen(result, "bac1")
    assert chosen["taxonomy_source"] == "silva"
    assert chosen["taxon_long"] == "Root;Bacteria;Pseudomonadota;Genus"
    assert chosen["phylum"] == "SilvaPhylum"
    assert chosen["pr2_taxon_long"] == "Root;Bacteria;Pr2;Genus"
    assert chosen["silva_domain"] == "Bacteria"
    assert not chosen["domain_disagreement"]


def test_agreed_eukaryota_keeps_pr2_organelle_suffix():
    pr2_taxon = "Root;Eukaryota:mito;Obazoa;Ophiostoma:mito;Ophiostoma_longirostellatum:mito"
    result = choose_taxonomy(
        _frame([_pr2("euk1", "Eukaryota:mito", 60, pr2_taxon, genus="Ophiostoma:mito")]),
        _frame([_silva("euk1", "Eukaryota", 99, "Root;Eukaryota;Fungi;Genus")]),
    )

    chosen = _chosen(result, "euk1")
    assert chosen["taxonomy_source"] == "pr2"
    assert chosen["taxon_long"] == pr2_taxon
    assert chosen["domain"] == "Eukaryota:mito"
    assert chosen["supergroup"] == "Pr2Supergroup"
    assert chosen["silva_taxon_long"] == "Root;Eukaryota;Fungi;Genus"
    assert not chosen["domain_disagreement"]


def test_agreed_archaea_keeps_silva():
    result = choose_taxonomy(
        _frame([_pr2("arc1", "Archaea", 80, "Root;Archaea;Pr2")]),
        _frame([_silva("arc1", "Archaea", 55, "Root;Archaea;Thermoproteota")]),
    )

    chosen = _chosen(result, "arc1")
    assert chosen["taxonomy_source"] == "silva"
    assert chosen["taxon_long"] == "Root;Archaea;Thermoproteota"
    assert not chosen["domain_disagreement"]


def test_disagreement_keeps_higher_domain_confidence():
    result = choose_taxonomy(
        _frame([
            _pr2("q_silva", "Eukaryota", 50, "Root;Eukaryota;Pr2"),
            _pr2("q_pr2", "Eukaryota:mito", 95, "Root;Eukaryota:mito;Pr2"),
        ]),
        _frame([
            _silva("q_silva", "Bacteria", 80, "Root;Bacteria;Silva"),
            _silva("q_pr2", "Bacteria", 20, "Root;Bacteria;SilvaLow"),
        ]),
    )

    silva_wins = _chosen(result, "q_silva")
    assert silva_wins["taxonomy_source"] == "silva"
    assert silva_wins["taxon_long"] == "Root;Bacteria;Silva"
    assert silva_wins["domain_disagreement"]
    assert silva_wins["pr2_taxon_long"] == "Root;Eukaryota;Pr2"

    pr2_wins = _chosen(result, "q_pr2")
    assert pr2_wins["taxonomy_source"] == "pr2"
    assert pr2_wins["taxon_long"] == "Root;Eukaryota:mito;Pr2"
    assert pr2_wins["domain_disagreement"]
    assert pr2_wins["silva_taxon_long"] == "Root;Bacteria;SilvaLow"


def test_tie_prefers_pr2_for_eukaryota_and_silva_for_bacteria():
    result = choose_taxonomy(
        _frame([
            _pr2("tie_euk", "Eukaryota", 70, "Root;Eukaryota;Pr2"),
            _pr2("tie_bac", "Unclassified", 50, "Root;Unclassified;Pr2"),
            _pr2("tie_other", "Foo", 10, "Root;Foo;Pr2"),
        ]),
        _frame([
            _silva("tie_euk", "Bacteria", 70, "Root;Bacteria;Silva"),
            _silva("tie_bac", "Bacteria", 50, "Root;Bacteria;Silva"),
            _silva("tie_other", "Bar", 10, "Root;Bar;Silva"),
        ]),
    )

    assert _chosen(result, "tie_euk")["taxonomy_source"] == "pr2"
    assert _chosen(result, "tie_bac")["taxonomy_source"] == "silva"
    assert _chosen(result, "tie_other")["taxonomy_source"] == "pr2"


def test_agreed_other_label_uses_higher_confidence_and_pr2_on_a_tie():
    result = choose_taxonomy(
        _frame([
            _pr2("higher_silva", "Unclassified", 10, "Root;Unclassified;Pr2"),
            _pr2("tie", "Unclassified", 40, "Root;Unclassified;Pr2Tie"),
        ]),
        _frame([
            _silva("higher_silva", "Unclassified", 80, "Root;Unclassified;Silva"),
            _silva("tie", "Unclassified", 40, "Root;Unclassified;SilvaTie"),
        ]),
    )

    assert _chosen(result, "higher_silva")["taxonomy_source"] == "silva"
    assert _chosen(result, "tie")["taxonomy_source"] == "pr2"
    assert not _chosen(result, "tie")["domain_disagreement"]


def test_query_present_in_only_one_database():
    result = choose_taxonomy(
        _frame([_pr2("only_pr2", "Eukaryota", 80, "Root;Eukaryota;Pr2")]),
        _frame([_silva("only_silva", "Bacteria", 80, "Root;Bacteria;Silva")]),
    )

    pr2_only = _chosen(result, "only_pr2")
    assert pr2_only["taxonomy_source"] == "pr2"
    assert not pr2_only["domain_disagreement"]
    assert pd.isna(pr2_only["silva_taxon_long"])
    assert pd.isna(pr2_only["silva_domain_confidence"])

    silva_only = _chosen(result, "only_silva")
    assert silva_only["taxonomy_source"] == "silva"
    assert pd.isna(silva_only["pr2_taxon_long"])
    assert pd.isna(silva_only["pr2_domain"])


def test_both_databases_empty():
    result = choose_taxonomy(pd.DataFrame(), pd.DataFrame())

    assert result.empty
    for column in (
        "query",
        "taxon_long",
        "confidence_long",
        "domain",
        "taxonomy_source",
        "domain_disagreement",
        "pr2_taxon_long",
        "pr2_domain",
        "pr2_domain_confidence",
        "silva_taxon_long",
        "silva_domain",
        "silva_domain_confidence",
    ):
        assert column in result.columns


def test_domain_confidence_falls_back_to_confidence_long():
    pr2 = pd.DataFrame([{
        "query": "from_long",
        "domain": "Eukaryota",
        "taxon_long": "Root;Eukaryota;Pr2",
        "confidence_long": "100;90;10",
    }])
    silva = pd.DataFrame([{
        "query": "from_long",
        "domain": "Bacteria",
        "taxon_long": "Root;Bacteria;Silva",
        "confidence_long": "100;15;10",
    }])

    chosen = _chosen(choose_taxonomy(pr2, silva), "from_long")
    assert chosen["taxonomy_source"] == "pr2"
    assert chosen["pr2_domain_confidence"] == 90
    assert chosen["silva_domain_confidence"] == 15


def test_abundance_domains_add_archaea_only_when_present():
    bacteria_and_euk = pd.DataFrame({
        "domain_clean": ["Bacteria", "Eukaryota"],
    })
    with_archaea = pd.DataFrame({
        "domain_clean": ["Bacteria", "Archaea"],
    })

    assert abundance_domains(bacteria_and_euk, "domain_clean") == ["Bacteria", "Eukaryota"]
    assert abundance_domains(with_archaea, "domain_clean") == ["Bacteria", "Eukaryota", "Archaea"]
