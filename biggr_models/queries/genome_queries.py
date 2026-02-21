from biggr_models.handlers import gene_handlers
from biggr_models.queries import utils
from cobradb.util import ref_tuple_to_str, ref_str_to_tuple
from cobradb.models import Genome, Chromosome, Model, GenomeRegion
from sqlalchemy import func, select

from sqlalchemy import inspect


def get_genomes_count(session, **kwargs):
    """Return the number of models in the database."""
    query = session.scalars(select(func.count(Genome.id))).first()
    return query


def get_all_genomes(session):
    """Get all genomes."""
    rows = session.execute(
        select(Genome.accession_type, Genome.accession_value).distinct()
    ).all()

    return [
        {"strain": r.accession_value, "accession_type": r.accession_type}
        for r in rows
        if r.accession_value is not None
    ]


def get_genomes(
    session,
    page=None,
    size=None,
    sort_column=None,
    sort_direction="ascending",
):
    # get the sort column
    columns = {
        "name": func.lower(Genome.accession_value),
        "organism": func.lower(Model.organism),
        "genome_type": func.lower(Genome.accession_type),
    }

    if sort_column is None:
        sort_column_object = None
    else:
        try:
            sort_column_object = columns[sort_column]
        except KeyError:
            print("Bad sort_column name: %s" % sort_column)
            sort_column_object = next(iter(columns.values()))

    query = select(Genome)
    query = utils._apply_order_limit_offset(
        query, sort_column_object, sort_direction, page, size
    )
    query = session.scalars(query).all()

    return [
        {
            "name": x.accession_value,
            "genome_ref_string": ref_tuple_to_str(x.accession_type, x.accession_value),
            "genome_type": x.accession_type,
            "organism": x.organism,
        }
        for x in query
    ]


def get_genome_and_models(genome_ref_string, session):
    accession_type, accession_value = ref_str_to_tuple(genome_ref_string)
    genome_db = session.scalars(
        select(Genome)
        .filter(Genome.accession_type == accession_type)
        .filter(Genome.accession_value == accession_value)
        .limit(1)
    ).first()
    models_db = session.scalars(
        select(Model).filter(Model.genome_id == genome_db.id)
    ).all()
    chromosomes_db = session.scalars(
        select(Chromosome).filter(Chromosome.genome_id == genome_db.id)
    ).all()

    genes_in_genome_url = f"/genomes/{genome_ref_string}/genes"
    genes_in_genome_columns = gene_handlers.GenesInGenomeListViewHandler.column_specs

    return {
        "name": genome_db.accession_value,
        "genome": genome_db,
        "genome_ref_string": ref_tuple_to_str(
            genome_db.accession_type, genome_db.accession_value
        ),
        "organism": genome_db.organism,
        "models": [x.bigg_id for x in models_db],
        "chromosomes": [x.ncbi_accession for x in chromosomes_db],
        "genes_in_genome_url": genes_in_genome_url,
        "genes_in_genome_columns": genes_in_genome_columns,
    }


def get_reactions_for_genome(genome_id, session):
    """Get all reactions associated with a genome through its models."""
    from cobradb.models import Model, ModelReaction, Reaction, UniversalReaction

    # Get reactions through the model associated with this genome
    reactions = session.execute(
        select(
            UniversalReaction.bigg_id,
            UniversalReaction.name,
            ModelReaction.gene_reaction_rule,
            ModelReaction.lower_bound,
            ModelReaction.upper_bound,
            ModelReaction.copy_number,
            ModelReaction.subsystem,
        )
        .join(Reaction, Reaction.universal_reaction_id == UniversalReaction.id)
        .join(ModelReaction, ModelReaction.reaction_id == Reaction.id)
        .join(Model, Model.id == ModelReaction.model_id)
        .filter(Model.genome_id == genome_id)
        .distinct()
    ).all()

    return [
        {
            "bigg_id": f"{r[0]}:{r[5]}" if r[5] != 1 else r[0],
            "name": r[1],
            "gene_reaction_rule": r[2],
            "lower_bound": r[3],
            "upper_bound": r[4],
            "copy_number": r[5],
            "subsystem": r[6],
        }
        for r in reactions
    ]


def get_metabolites_for_genome(genome_id, session):
    """Get all metabolites associated with a genome through its models."""
    from cobradb.models import (
        Model,
        ModelCompartmentalizedComponent,
        CompartmentalizedComponent,
        Component,
    )

    # Get metabolites through the model associated with this genome
    metabolites = session.execute(
        select(
            Component.bigg_id,
            Component.name,
            Component.formula,
            Component.charge,
            CompartmentalizedComponent.bigg_id.label("compartmentalized_bigg_id"),
        )
        .join(
            CompartmentalizedComponent,
            CompartmentalizedComponent.component_id == Component.id,
        )
        .join(
            ModelCompartmentalizedComponent,
            ModelCompartmentalizedComponent.compartmentalized_component_id
            == CompartmentalizedComponent.id,
        )
        .join(Model, Model.id == ModelCompartmentalizedComponent.model_id)
        .filter(Model.genome_id == genome_id)
        .distinct()
    ).all()

    return [
        {
            "bigg_id": m[0],
            "name": m[1],
            "formula": m[2],
            "charge": m[3],
            "compartmentalized_bigg_id": m[4],
        }
        for m in metabolites
    ]


def get_genomes_with_chromosomes(
    accession_id,
    session,
    gene_id_filter=None,
    include_metabolites=True,
    include_reactions=True,
):
    if not accession_id:
        return []

    genomes = session.scalars(
        select(Genome).filter(Genome.accession_value == accession_id)
    ).all()

    if not genomes:
        return []

    genome_ids = [g.id for g in genomes]

    chromosomes = session.scalars(
        select(Chromosome).filter(Chromosome.genome_id.in_(genome_ids))
    ).all()
    if not chromosomes:
        chrom_map = {}
    else:
        chrom_ids = [c.id for c in chromosomes]

        region_query = select(GenomeRegion).filter(
            GenomeRegion.chromosome_id.in_(chrom_ids)
        )

        if gene_id_filter is not None:
            filt_ids = set(map(int, gene_id_filter))
            if filt_ids:
                region_query = region_query.filter(GenomeRegion.id.in_(filt_ids))
            else:
                region_query = region_query.filter(False)

        regions = session.scalars(region_query).all()

        region_map = {}
        for r in regions:
            region_dict = {
                col.key: getattr(r, col.key)
                for col in inspect(GenomeRegion).mapper.column_attrs
            }
            region_map.setdefault(r.chromosome_id, []).append(region_dict)

        chrom_map = {}
        for c in chromosomes:
            chrom_dict = {
                col.key: getattr(c, col.key)
                for col in inspect(Chromosome).mapper.column_attrs
            }
            chrom_dict["genome_region"] = region_map.get(c.id, [])
            chrom_map.setdefault(c.genome_id, []).append(chrom_dict)

    results = []
    for g in genomes:
        genome_dict = {
            col.key: getattr(g, col.key) for col in inspect(Genome).mapper.column_attrs
        }
        genome_dict["chromosome"] = chrom_map.get(g.id, [])

        # Add metabolites if requested
        if include_metabolites:
            genome_dict["metabolites"] = get_metabolites_for_genome(g.id, session)

        # Add reactions if requested
        if include_reactions:
            genome_dict["reactions"] = get_reactions_for_genome(g.id, session)

        results.append(genome_dict)

    return results
