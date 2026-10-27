import os
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt


ORGANELLE_SUFFIXES = {
    "nucl": "nucleomorph",
    "plas": "plastid",
    "apic": "apicoplast",
    "chrom": "chromatophore",
    "mito": "mitochondrion",
}

ORIGIN_COLORS = {
    "cellular": "grey",
    "nuclear": "skyblue",
    "mito": "orange",
    "plas": "green",
    "nucl": "purple",
    "apic": "red",
    "chrom": "brown",
}

def extract_genus_and_type(taxon_long: str):
    """
    Extract genus and origin type from PR2 taxonomy string.

    Example:
    Root;Eukaryota:mito;...;Ophiostoma:mito;Ophiostoma_longirostellatum:mito
    -> genus='Ophiostoma', origin_type='mito'
    """

    parts = taxon_long.split(";")
    if len(parts) < 2:
        return ("Unknown", "nuclear")

    genus_part = parts[-2]

    if ":" in genus_part:
        genus, suffix = genus_part.split(":", 1)

        # validate suffix
        if suffix in ORGANELLE_SUFFIXES:
            origin_type = suffix
        else:
            origin_type = "nuclear"

    else:
        genus = genus_part
        origin_type = "nuclear"

    return genus, origin_type


def cleaned_domain(series: pd.Series) -> pd.Series:
    """Drop a PR2 organelle suffix so Eukaryota:mito compares as Eukaryota."""
    return series.fillna("").astype(str).str.replace(r":.*", "", regex=True)


def abundance_domains(df: pd.DataFrame, domain_column: str) -> list[str]:
    """Bacteria and Eukaryota always. Archaea only when a chosen row is archaeal."""
    domains = ["Bacteria", "Eukaryota"]
    if domain_column not in df.columns or df.empty:
        return domains
    if (cleaned_domain(df[domain_column]) == "Archaea").any():
        domains.append("Archaea")
    return domains


def _is_missing(value) -> bool:
    if value is None:
        return True
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        return False
    return bool(missing)


def _domain_label(row):
    """Domain string as reported, including a PR2 organelle suffix when present."""
    if row is None:
        return pd.NA
    if "domain" in row.index and not _is_missing(row["domain"]):
        return row["domain"]
    taxon = row["taxon_long"] if "taxon_long" in row.index else None
    if isinstance(taxon, str) and taxon.count(";") >= 1:
        return taxon.split(";")[1]
    return pd.NA


def _cleaned_domain_label(row):
    label = _domain_label(row)
    if _is_missing(label):
        return None
    text = str(label).strip()
    if not text or text.lower() == "nan":
        return None
    return text.split(":", 1)[0]


def _domain_confidence(row) -> float:
    """Numeric domain confidence. Missing scores lose to any real score."""
    if row is None:
        return float("-inf")
    if "domain_confidence" in row.index and not _is_missing(row["domain_confidence"]):
        try:
            return float(row["domain_confidence"])
        except (TypeError, ValueError):
            pass
    confidence_long = row["confidence_long"] if "confidence_long" in row.index else None
    if isinstance(confidence_long, str):
        parts = confidence_long.split(";")
        if len(parts) > 1 and parts[1] != "":
            try:
                return float(parts[1])
            except ValueError:
                return float("-inf")
    return float("-inf")


def _recorded_confidence(row):
    if row is None:
        return pd.NA
    score = _domain_confidence(row)
    if score == float("-inf"):
        return pd.NA
    return score


def _taxon_long(row):
    if row is None or "taxon_long" not in row.index or _is_missing(row["taxon_long"]):
        return pd.NA
    return row["taxon_long"]


def _index_by_query(df: pd.DataFrame) -> dict:
    if df is None or df.empty or "query" not in df.columns:
        return {}
    indexed = {}
    for _, row in df.iterrows():
        if _is_missing(row["query"]):
            continue
        query = str(row["query"])
        current = indexed.get(query)
        # A repeated query keeps the stronger domain call.
        if current is None or _domain_confidence(row) > _domain_confidence(current):
            indexed[query] = row
    return indexed


def _choose_source(pr2_row, silva_row) -> str:
    """SILVA for Bacteria and Archaea, PR2 for Eukaryota. Confidence breaks disagreements."""
    if pr2_row is None:
        return "silva"
    if silva_row is None:
        return "pr2"

    pr2_domain = _cleaned_domain_label(pr2_row)
    silva_domain = _cleaned_domain_label(silva_row)
    pr2_conf = _domain_confidence(pr2_row)
    silva_conf = _domain_confidence(silva_row)

    if pr2_domain == silva_domain:
        if pr2_domain in {"Bacteria", "Archaea"}:
            return "silva"
        if pr2_domain == "Eukaryota":
            return "pr2"
        if silva_conf > pr2_conf:
            return "silva"
        return "pr2"

    if pr2_conf > silva_conf:
        return "pr2"
    if silva_conf > pr2_conf:
        return "silva"
    # Equal confidence: keep the database that owns the domain it called.
    if pr2_domain == "Eukaryota":
        return "pr2"
    if silva_domain in {"Bacteria", "Archaea"}:
        return "silva"
    return "pr2"


def _empty_taxonomy() -> pd.DataFrame:
    return pd.DataFrame(columns=[
        "query",
        "confidence_long",
        "taxon_long",
        "domain",
        "taxonomy_source",
        "domain_disagreement",
        "pr2_taxon_long",
        "pr2_domain",
        "pr2_domain_confidence",
        "silva_taxon_long",
        "silva_domain",
        "silva_domain_confidence",
    ])


def _domains_disagree(pr2_row, silva_row) -> bool:
    if pr2_row is None or silva_row is None:
        return False
    return _cleaned_domain_label(pr2_row) != _cleaned_domain_label(silva_row)


def choose_taxonomy(pr2_df: pd.DataFrame, silva_df: pd.DataFrame) -> pd.DataFrame:
    """One row per query, with the other database's call kept beside it."""
    pr2_by_query = _index_by_query(pr2_df)
    silva_by_query = _index_by_query(silva_df)
    queries = list(pr2_by_query) + [query for query in silva_by_query if query not in pr2_by_query]
    if not queries:
        return _empty_taxonomy()

    records = []
    for query in queries:
        pr2_row = pr2_by_query.get(query)
        silva_row = silva_by_query.get(query)
        source = _choose_source(pr2_row, silva_row)
        chosen = pr2_row if source == "pr2" else silva_row
        record = chosen.to_dict()
        if _is_missing(record.get("domain")):
            label = _domain_label(chosen)
            if not _is_missing(label):
                record["domain"] = label
        record["taxonomy_source"] = source
        record["domain_disagreement"] = _domains_disagree(pr2_row, silva_row)
        record["pr2_taxon_long"] = _taxon_long(pr2_row)
        record["pr2_domain"] = _domain_label(pr2_row)
        record["pr2_domain_confidence"] = _recorded_confidence(pr2_row)
        record["silva_taxon_long"] = _taxon_long(silva_row)
        record["silva_domain"] = _domain_label(silva_row)
        record["silva_domain_confidence"] = _recorded_confidence(silva_row)
        records.append(record)
    return pd.DataFrame.from_records(records)


def read_decipher_table(path: str) -> pd.DataFrame:
    if os.path.exists(path):
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame(columns=["query", "confidence_long", "taxon_long", "domain"])


def plot_genus_abundance(df: pd.DataFrame, out_pdf: Path, top_n=20):
    """
    Plot total read abundances per genus, stacked by origin_type.
    Facet by domain (Bacteria / Eukaryota, plus Archaea when present).
    """

    domains = abundance_domains(df, "domain")
    cleaned = cleaned_domain(df["domain"])

    fig, axes = plt.subplots(
        nrows=len(domains),
        ncols=1,
        figsize=(10, 3 * len(domains)),
        sharex=True)

    if len(domains) == 1:
        axes = [axes]

    for ax, domain in zip(axes, domains):

        df_dom = df[cleaned == domain]

        if df_dom.empty:
            ax.set_title(domain)
            ax.axis("off")
            continue

        counts = (
            df_dom.groupby(["genus", "origin_type"])
            .size()
            .unstack(fill_value=0))

        total_counts = counts.sum(axis=1)

        if len(total_counts) > top_n:
            top_genera = total_counts.nlargest(top_n).index
            counts_top = counts.loc[top_genera]
            counts_other = counts.drop(top_genera).sum()
            counts_top.loc["Other"] = counts_other
            counts = counts_top

        counts = counts.loc[counts.sum(axis=1).sort_values().index]

        counts.plot(
            kind="barh",
            stacked=True,
            ax=ax,
            color=ORIGIN_COLORS,
            edgecolor="black")

        ax.set_title(domain)
        ax.set_ylabel("Genus")

    axes[-1].set_xlabel("Number of reads")

    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()

def plot_taxon_abundance(df, rank, out_pdf, top_n=20):

    domains = abundance_domains(df, "domain_clean")

    fig, axes = plt.subplots(
        nrows=len(domains),
        ncols=1,
        figsize=(10, 3 * len(domains)),
        sharex=True
    )

    if len(domains) == 1:
        axes = [axes]

    for ax, domain in zip(axes, domains):

        df_dom = df[df["domain_clean"] == domain]

        if df_dom.empty:
            ax.axis("off")
            continue

        counts = (
            df_dom.groupby([rank, "origin_type"])
            .size()
            .unstack(fill_value=0)
        )

        totals = counts.sum(axis=1)

        if len(totals) > top_n:
            top = totals.nlargest(top_n).index
            counts_top = counts.loc[top]
            counts_other = counts.drop(top).sum()
            counts_top.loc["Other"] = counts_other
            counts = counts_top

        counts = counts.loc[counts.sum(axis=1).sort_values().index]

        counts.plot(
            kind="barh",
            stacked=True,
            ax=ax,
            color=ORIGIN_COLORS,
            edgecolor="black"
        )

        ax.set_title(domain)
        ax.set_ylabel(rank.capitalize())

    axes[-1].set_xlabel("Number of reads")

    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()

def run(ctx):
    mapping_file = ctx.require("mapping")
    out_file = ctx.artifact("taxa_ssu_reads", "taxonomy", ".tsv")

    mapping = pd.read_csv(mapping_file, sep="\t")
    pr2_df = read_decipher_table(os.path.join(ctx.output_dir, "decipher_pr2.tsv"))
    silva_df = read_decipher_table(os.path.join(ctx.output_dir, "decipher_silva.tsv"))
    df = choose_taxonomy(pr2_df, silva_df)

    # normalize origin to merge with mapping RNAME
    df['origin'] = (
        df['query']
        .str.replace(r'_SSU.*','',regex=True)
        .str.split().str[0]
        .str.replace(r':\d+-\d+\([+-]\)\|.*$','',regex=True))

    contig_tax = (
        df.sort_values(['origin','confidence_long'], ascending=[True, False])
          .groupby('origin', as_index=False)
          .first())

    df_annotated_reads = (
        mapping[['QNAME','RNAME','FLAG','POS','MAPQ']]
        .merge(contig_tax, left_on='RNAME', right_on='origin', how='left')
        .dropna(subset=['taxon_long'])
    )

    # extract genus and origin type
    df_annotated_reads[['genus', 'origin_type']] = df_annotated_reads['taxon_long'].apply(
        lambda x: pd.Series(extract_genus_and_type(x))
    )

    df_annotated_reads["domain_clean"] = df_annotated_reads["domain"].str.replace(r":.*", "", regex=True)

    df_annotated_reads.loc[
        df_annotated_reads["domain_clean"] != "Eukaryota",
        "origin_type"
    ] = "cellular" 


    df_annotated_reads.to_csv(out_file, sep="\t", index=False)
    ctx.register("taxa_ssu_reads", out_file)

    # generate abundance plot
    pdf_out = Path(mapping_file).parent / "genus_abundance.pdf"
    plot_genus_abundance(df_annotated_reads, pdf_out)
    ctx.log(f"Genus abundance plot written to {pdf_out}")
    ctx.register("genus_abundance_plot", pdf_out)

    tax_ranks = ["class", "order", "family", "genus"]

    plot_dir = Path(mapping_file).parent

    for rank in tax_ranks:

        pdf_out = plot_dir / f"{rank}_abundance.pdf"

        plot_taxon_abundance(
            df_annotated_reads,
            rank=rank,
            out_pdf=pdf_out
        )

    ctx.log(f"{rank} abundance plot written to {pdf_out}")