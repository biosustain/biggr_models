from typing import Dict, Any
from biggr_models.queries import utils
from cobradb.util import ref_tuple_to_str
from cobradb.models import (
    Chromosome,
    Gene,
    GeneReactionMatrix,
    Genome,
    Model,
    ModelGene,
    ModelReaction,
    Reaction,
    UniversalReaction,
    GenomeRegion,
)

from sqlalchemy import func, select
from sqlalchemy.orm import Session, subqueryload, contains_eager


def get_gene_ids_for_gene_name(name, session):
    """Get the gene ids for a gene name."""
    return session.scalars(
        select(Gene.id).filter(func.lower(Gene.name) == func.lower(name))
    ).all()


def get_genes(gene_ids, session):
    """Get the genes for a list of gene ids."""
    rows = session.execute(
        select(
            Gene.id, Gene.bigg_id, Gene.name, Gene.locus_tag, Gene.mapped_to_genbank
        ).filter(Gene.id.in_(list(gene_ids)))
    ).all()

    return [dict(r._mapping) for r in rows]


def get_all_genes(session):
    """Get all genes."""
    rows = session.execute(
        select(Gene.id, Gene.bigg_id, Gene.name, Gene.locus_tag)
    ).all()

    return [dict(r._mapping) for r in rows]


def get_all_genes_with_urls(session):
    """Get all genes with their associated genome and model URLs."""
    # Genome URL: /genomes/{accession_type}:{accession_value}/genes/{bigg_id}
    genome_rows = session.execute(
        select(Gene.name, Gene.bigg_id, Genome.accession_type, Genome.accession_value)
        .join(Chromosome, Chromosome.id == Gene.chromosome_id)
        .join(Genome, Genome.id == Chromosome.genome_id)
        .filter(Gene.name.isnot(None))
    ).all()

    # Model URL: /models/{model_bigg_id}/genes/{bigg_id}
    model_rows = session.execute(
        select(Gene.name, Gene.bigg_id, Model.bigg_id)
        .join(ModelGene, ModelGene.gene_id == Gene.id)
        .join(Model, Model.id == ModelGene.model_id)
        .filter(Gene.name.isnot(None))
    ).all()

    # Group URLs by gene name
    gene_map = {}
    for name, bigg_id, acc_type, acc_value in genome_rows:
        gene_map.setdefault(name, set()).add(
            f"/genomes/{acc_type}:{acc_value}/genes/{bigg_id}"
        )
    for name, bigg_id, model_bigg_id in model_rows:
        gene_map.setdefault(name, set()).add(
            f"/models/{model_bigg_id}/genes/{bigg_id}"
        )

    return [
        {"name": name, "urls": sorted(urls)}
        for name, urls in sorted(gene_map.items())
    ]


def get_all_gene_strain_pairs(session):
    """Get all gene-strain pairs with their genome URLs."""
    rows = session.execute(
        select(
            Gene.name,
            Gene.bigg_id,
            Gene.locus_tag,
            Genome.accession_type,
            Genome.accession_value,
        )
        .join(Chromosome, Chromosome.id == Gene.chromosome_id)
        .join(Genome, Genome.id == Chromosome.genome_id)
        .filter(Gene.name.isnot(None))
    ).all()

    # Group by (gene_name, strain), collecting URLs
    pair_map = {}
    for name, bigg_id, locus_tag, acc_type, acc_value in rows:
        key = (name, acc_value)
        if key not in pair_map:
            pair_map[key] = {"gene": name, "strain": acc_value, "urls": []}
        pair_map[key]["urls"].append(
            f"/genomes/{acc_type}:{acc_value}/genes/{bigg_id}"
        )

    return list(pair_map.values())


def get_genes_with_genome_region(gene_ids, session):
    """Get gene + genome_region flattened data with genome_gene_url for a list of gene IDs."""
    gene_ids = list(gene_ids)
    rows = session.execute(
        select(
            Gene.id,
            Gene.name,
            Gene.locus_tag,
            Gene.mapped_to_genbank,
            GenomeRegion.bigg_id,
            GenomeRegion.chromosome_id,
            GenomeRegion.leftpos,
            GenomeRegion.rightpos,
            GenomeRegion.strand,
            GenomeRegion.type,
            GenomeRegion.dna_sequence,
            GenomeRegion.protein_sequence,
            Genome.accession_type,
            Genome.accession_value,
        )
        .join(GenomeRegion, GenomeRegion.id == Gene.id)
        .join(Chromosome, Chromosome.id == GenomeRegion.chromosome_id)
        .join(Genome, Genome.id == Chromosome.genome_id)
        .filter(Gene.id.in_(gene_ids))
    ).all()

    return [
        {
            "id": r[0],
            "name": r[1],
            "locus_tag": r[2],
            "mapped_to_genbank": r[3],
            "bigg_id": r[4],
            "chromosome_id": r[5],
            "leftpos": r[6],
            "rightpos": r[7],
            "strand": r[8],
            "type": r[9],
            "dna_sequence": r[10],
            "protein_sequence": r[11],
            "genome_gene_url": f"/genomes/{r[12]}:{r[13]}/genes/{r[4]}",
        }
        for r in rows
    ]


def get_model_genes_for_gene_ids(gene_ids, session):
    """Get model_gene + model flattened data with model_gene_url for a list of gene IDs."""
    gene_ids = list(gene_ids)
    rows = session.execute(
        select(
            ModelGene.id,
            ModelGene.gene_id,
            ModelGene.model_id,
            Model.bigg_id,
            Model.organism,
            Model.genome_id,
            Model.published_filename,
            Gene.bigg_id,
        )
        .join(Model, Model.id == ModelGene.model_id)
        .join(Gene, Gene.id == ModelGene.gene_id)
        .filter(ModelGene.gene_id.in_(gene_ids))
    ).all()

    return [
        {
            "model_gene_id": r[0],
            "gene_id": r[1],
            "model_id": r[2],
            "bigg_id": r[3],
            "organism": r[4],
            "genome_id": r[5],
            "published_filename": r[6],
            "model_gene_url": f"/models/{r[3]}/genes/{r[7]}",
        }
        for r in rows
    ]


def get_urls_for_gene_ids(gene_ids, session):
    """Get genome and model URLs for a list of gene IDs, grouped by gene ID."""
    gene_ids = list(gene_ids)

    # Genome URL: /genomes/{accession_type}:{accession_value}/genes/{bigg_id}
    genome_rows = session.execute(
        select(Gene.id, Gene.bigg_id, Genome.accession_type, Genome.accession_value)
        .join(Chromosome, Chromosome.id == Gene.chromosome_id)
        .join(Genome, Genome.id == Chromosome.genome_id)
        .filter(Gene.id.in_(gene_ids))
    ).all()

    # Model URL: /models/{model_bigg_id}/genes/{bigg_id}
    model_rows = session.execute(
        select(Gene.id, Gene.bigg_id, Model.bigg_id)
        .join(ModelGene, ModelGene.gene_id == Gene.id)
        .join(Model, Model.id == ModelGene.model_id)
        .filter(Gene.id.in_(gene_ids))
    ).all()

    url_map = {}
    for gid, bigg_id, acc_type, acc_value in genome_rows:
        url_map.setdefault(gid, set()).add(
            f"/genomes/{acc_type}:{acc_value}/genes/{bigg_id}"
        )
    for gid, bigg_id, model_bigg_id in model_rows:
        url_map.setdefault(gid, set()).add(
            f"/models/{model_bigg_id}/genes/{bigg_id}"
        )

    return {gid: sorted(urls) for gid, urls in url_map.items()}


def get_genome_region_for_gene_id(ids, session):
    """Get the genome region for a gene id."""
    rows = session.execute(
        select(
            GenomeRegion.id,
            GenomeRegion.chromosome_id,
            GenomeRegion.bigg_id,
            GenomeRegion.leftpos,
            GenomeRegion.rightpos,
            GenomeRegion.strand,
            GenomeRegion.type,
            GenomeRegion.dna_sequence,
            GenomeRegion.protein_sequence,
        ).filter(GenomeRegion.id.in_(list(ids)))
    ).all()

    return [dict(r._mapping) for r in rows]


def get_model_genes_count(model_bigg_id, session):
    """Get the number of gene for the given model."""
    return session.scalars(
        session.query(func.count(Gene.id))
        .join(ModelGene)
        .join(Model)
        .filter(Model.bigg_id == model_bigg_id)
    ).first()


def get_model_genes(
    model_bigg_id,
    session,
    page=None,
    size=None,
    sort_column=None,
    sort_direction="ascending",
    **kwargs
):
    """Get model genes.

    Arguments
    ---------

    model_bigg_id: The bigg id of the model to retrieve genes.

    session: An ome session object.

    page: The page, or None for all pages.

    size: The page length, or None for all pages.

    sort_column: The name of the column to sort. Must be one of 'bigg_id', 'name',
    'model_bigg_id', and 'organism'.

    sort_direction: Either 'ascending' or 'descending'.

    Returns
    -------

    A list of objects with keys 'bigg_id', 'name', 'model_bigg_id', and
    'organism'.

    """
    # get the sort column
    columns = {
        "bigg_id": func.lower(Gene.bigg_id),
        "name": func.lower(Gene.name),
        "model_bigg_id": func.lower(Model.bigg_id),
        "organism": func.lower(Model.organism),
    }

    if sort_column is None:
        sort_column_object = None
    else:
        try:
            sort_column_object = columns[sort_column]
        except KeyError:
            raise ValueError("Bad sort_column name: %s" % sort_column)

    # set up the query
    query = (
        select(Gene.bigg_id, Gene.name, Model.bigg_id, Model.organism)
        .join(ModelGene, ModelGene.gene_id == Gene.id)
        .join(Model, Model.id == ModelGene.model_id)
        .filter(Model.bigg_id == model_bigg_id)
    )

    # order and limit
    query = utils._apply_order_limit_offset(
        query, sort_column_object, sort_direction, page, size
    )

    query = session.execute(query).all()
    return [
        {"bigg_id": x[0], "name": x[1], "model_bigg_id": x[2], "organism": x[3]}
        for x in query
    ]


def get_model_gene(gene_bigg_id, model_bigg_id, session):
    result_db = session.execute(
        select(
            Gene.bigg_id,
            Gene.name,
            Gene.leftpos,
            Gene.rightpos,
            Model.bigg_id,
            Gene.strand,
            Chromosome.ncbi_accession,
            Genome.accession_type,
            Genome.accession_value,
            Gene.mapped_to_genbank,
            Gene.dna_sequence,
            Gene.protein_sequence,
        )
        .join(ModelGene, ModelGene.gene_id == Gene.id)
        .join(Model, Model.id == ModelGene.model_id)
        .outerjoin(Genome, Genome.id == Model.genome_id)
        .outerjoin(Chromosome, Chromosome.id == Gene.chromosome_id)
        .filter(Gene.bigg_id == gene_bigg_id)
        .filter(Model.bigg_id == model_bigg_id)
        .limit(1)
    ).first()
    if result_db is None:
        raise utils.NotFoundError(
            "Gene %s not found in model %s" % (gene_bigg_id, model_bigg_id)
        )

    reaction_db = session.execute(
        select(
            Reaction.bigg_id,
            ModelReaction.gene_reaction_rule,
            UniversalReaction.name,
            UniversalReaction.bigg_id,
        )
        .join(UniversalReaction, UniversalReaction.id == Reaction.universal_reaction_id)
        .join(ModelReaction, ModelReaction.reaction_id == Reaction.id)
        .join(Model, Model.id == ModelReaction.model_id)
        .join(
            GeneReactionMatrix, GeneReactionMatrix.model_reaction_id == ModelReaction.id
        )
        .join(ModelGene, ModelGene.id == GeneReactionMatrix.model_gene_id)
        .join(Gene, Gene.id == ModelGene.gene_id)
        .filter(Model.bigg_id == model_bigg_id)
        .filter(Gene.bigg_id == gene_bigg_id)
    )
    reaction_results = [
        {"bigg_id": r[3], "reaction_id": r[0], "gene_reaction_rule": r[1], "name": r[2]}
        for r in reaction_db
    ]

    # synonym_db = id_queries._get_db_links_for_model_gene(gene_bigg_id, session)
    synonym_db = {}

    # old_id_results = id_queries._get_old_ids_for_model_gene(
    #     gene_bigg_id, model_bigg_id, session
    # )
    old_id_results = []

    return {
        "bigg_id": result_db[0],
        "name": result_db[1],
        "leftpos": result_db[2],
        "rightpos": result_db[3],
        "model_bigg_id": result_db[4],
        "strand": result_db[5],
        "chromosome_ncbi_accession": result_db[6],
        "genome_ref_string": ref_tuple_to_str(result_db[7], result_db[8]),
        "genome_name": result_db[8],
        "mapped_to_genbank": result_db[9],
        "dna_sequence": result_db[10],
        "protein_sequence": result_db[11],
        "reactions": reaction_results,
        "database_links": synonym_db,
        "old_identifiers": old_id_results,
    }


def get_gene(
    session: Session, accession_type: str, accession_value: str, gene_bigg_id: str
) -> Dict[str, Any]:
    gene_db = session.scalars(
        select(Gene)
        .options(
            contains_eager(Gene.chromosome).contains_eager(Chromosome.genome),
            subqueryload(Gene.model_genes).joinedload(ModelGene.model),
        )
        .join(Gene.chromosome)
        .join(Chromosome.genome)
        .filter(Genome.accession_type == accession_type)
        .filter(Genome.accession_value == accession_value)
        .filter(Gene.bigg_id == gene_bigg_id)
    ).first()
    if gene_db is None:
        raise utils.NotFoundError("Could not find gene")

    return {"gene": gene_db}
