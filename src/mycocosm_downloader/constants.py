# TaxIds that need manual correction (merged or outdated in NCBI)
ALT_TAXID: dict[str, str] = {
    "1140396": "2735509",
    "585595": "2587404",
    "5145": "2587412",
    "1437435": "2211651",
}

# Main branches of the JGI fungal taxonomy tree
# See https://genome.jgi.doe.gov/programs/fungi/index.jsf
JGI_TREE_BRANCHES: set[str] = {
    "pucciniomycotina",
    "ustilaginomycotina",
    "agaricomycetes",
    "dacrymycetes",
    "tremellomycetes",
    "wallemiomycetes",
    "pezizomycetes",
    "orbiliomycetes",
    "eurotiomycetes",
    "dothideomycetes",
    "lecanoromycetes",
    "leotiomycetes",
    "sordariomycetes",
    "xylonomycetes",
    "saccharomycotina",
    "taphrinomycotina",
    "glomeromycotina",
    "mortierellomycotina",
    "mucoromycotina",
    "zoopagomycotina",
    "entomophthoromycotina",
    "kickxellomycotina",
    "blastocladiomycota",
    "chytridiomycetes",
    "monoblepharidomycetes",
    "neocallimastigomycota",
    "microsporidia",
    "cryptomycota",
}

# Metaprojects or old versions to exclude from downloads
PORTALS_TO_REMOVE: set[str] = {
    "Aciri1_meta",
    "Altbr1",
    "Pospl1",  # superseded by PosplRSB12_1
    "Rhoto_IFO0880_2",  
}

# Translates project shortname to organism name (encoding fixes for some species)
MANUAL_PROJECT_NAMES: dict[str, str] = {
    "Boledp1": "Boletus edulis Přilba v1.0",
    "Rusoch1": "Russula ochroleuca Přilba v1.0",
    "Amarub1": "Amanita rubescens Přilba v1.0",
    "Ruseme1": "Russula emetica Přilba v1.0",
}

# Assembly filenames to ignore (not real assemblies, old versions, or meta-samples)
EXCLUDE_ASSEMBLIES: set[str] = {
    "1034997.Tuber_borchii_Tbo3840.standard.main.scaffolds.fasta.gz",
    "Spofi1.draft.mito.scaffolds.fasta.gz",
    "Patat1.draft.mito.scaffolds.fasta.gz",
    "PleosPC9_1_Assembly_scaffolds.fasta.gz",
    "Neuhi1_PlasmidAssemblyScaffolds.fasta.gz",
    "CocheC5_1_assembly_scaffolds.fasta.gz",
    "Alternaria_brassicicola_masked_assembly.fasta.gz",
    "Aciri1_meta_AssemblyScaffolds.fasta.gz",
    "Rhoto_IFO0880_2_AssemblyScaffolds.fasta.gz",
    "StenotrophomonasSp_AssemblyScaffolds.fasta.gz",
    "PseudomonasSp_AssemblyScaffolds.fasta.gz",
    "EurotioJF034F_1_RiboAssemblyScaffolds.fasta.gz",
}

# GFF filenames to ignore (have "GeneCatalog" magic string but are unwanted)
IGNORE_GFF_FILENAMES: set[str] = {
    "Aciri1_meta_GeneCatalog_genes_20111216.gff.gz",
    "Exoaq1_GeneCatalog_20160901.gff3.gz",
    "Exoaq1_GeneCatalog_20160828.gff3.gz",
    "Fonpe1_GeneCatalog_20160901.gff3.gz",
    "Copmic2_FM1_removed_alleles.gff.gz",
}

# US timezone abbreviations to UTC offsets (for timestamp parsing)
US_TIMEZONES: dict[str, str] = {
    "EST": "-0500",
    "CST": "-0600",
    "MST": "-0700",
    "PST": "-0800",
    "EDT": "-0400",
    "CDT": "-0500",
    "MDT": "-0600",
    "PDT": "-0700",
}

# Weekday abbreviations to numeric values (for locale-independent parsing)
US_DAY_NUMBERS: dict[str, str] = {
    "Sun": "0",
    "Mon": "1",
    "Tue": "2",
    "Wed": "3",
    "Thu": "4",
    "Fri": "5",
    "Sat": "6",
}

# Month abbreviations to numeric values (for locale-independent parsing)
US_MONTH_NUMBERS: dict[str, str] = {
    "Jan": "01",
    "Feb": "02",
    "Mar": "03",
    "Apr": "04",
    "May": "05",
    "Jun": "06",
    "Jul": "07",
    "Aug": "08",
    "Sep": "09",
    "Oct": "10",
    "Nov": "11",
    "Dec": "12",
}