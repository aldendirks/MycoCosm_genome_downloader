#!/usr/bin/env python

"""
This script downloads all fungal genomes from JGI: 
genome assembly (fasta) + gene annotations (gff)

Also, get a summary of metadata:
TaxId   JGI_short_name  Complete name   Original name   +files location

Note that JGI only lists the latest version (i.e. if v2.0 is available for a 
certain assembly, v1.0 will not be available anymore)
"""

import argparse
import codecs
import csv
import shutil
import subprocess as sp
import time
from collections import defaultdict
from datetime import datetime
from getpass import getpass
from pathlib import Path
from subprocess import DEVNULL, STDOUT, CalledProcessError

from ete4.ncbi_taxonomy import NCBITaxa
from lxml import etree

__author__ = "Jorge Navarro, Alden Dirks"
__version__ = "v1.4.0"


# ============================================================================
# CONSTANTS
# ============================================================================

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


# ============================================================================
# DATA STRUCTURES
# ============================================================================

class JGI_Project:
    """
    Store information about each JGI sequencing project (portal).
    """
    def __init__(self):
        self.TaxId = ""  # NCBI TaxId (not lineage)
        self.lineage_set = set()  # to quickly check major fungal taxonomy branches
        self.lineage_list = []

        self.portal = ""  # Fam[:2] + Sp[:1] + mystery number. a.k.a.: "shortname"
        self.name = ""  # Portal full name. Include things like "v2.0"
        self.org_name = ""  # Portal name without strings like "v2.0"

        self.assembly_file = ""
        self.assembly_url = ""
        self.assembly_size = 0
        self.gff_file = ""
        self.gff_url = ""
        self.gff_timestamp = ""
        self.gff_size = 0

        self.project_path = ""  # Path object

        self.date = ""  # Ideally, a Time object

        self.is_restricted = False

        return


# ============================================================================
# COMMAND-LINE INTERFACE
# ============================================================================

def command_parser():
    parser = argparse.ArgumentParser(
        description=f"Downloads fungal genomes from JGI's MycoCosm.\
            Version: {__version__}."
    )

    group_get = parser.add_argument_group("get necessary data")
    group_get.add_argument(
        "-u",
        "--update",
        help="Updates NCBI Taxonomy database and exits. Used to infer \
            lineage from TaxId.",
        default=False,
        action="store_true",
    )
    group_get.add_argument(
        "--getgenomelist",
        help="Gets a fresh download of MycoCosm's genome list \
            'MycoCosm_Genome_list.csv' and exits. Used to get metadata \
            of each project (TaxId, full name). Stored in the folder \
            specified by --outputfolder.",
        default=False,
        action="store_true",
    )
    group_get.add_argument(
        "--getxml",
        help="Builds an XML file with all files per project. Stored in the \
            folder specified by --outputfolder. If already present, only \
            missing data will be added. Uses 'MycoCosm_Genome_list.csv' (which \
            should exist in the output folder).",
        default=False,
        action="store_true",
    )
    group_get.add_argument(
        "--getprevious",
        type=Path,
        help="Gets a list of locations of .gz files, downloaded from a \
            previous run. Argument to this parameter is the base folder to \
            start looking recursively for .gz files. File will be stored in \
            the same place as this script (name: \
            'previously_downloaded_files.tsv')."
    )

    group_input = parser.add_argument_group("input data")
    group_input.add_argument(
        "-c",
        "--csv",
        help="List of JGI genomes with some metadata like TaxId.",
        type=Path
    )
    group_input.add_argument(
        "-x", 
        "--xml", 
        help="JGI file list with file URL paths.",
        type=Path
    )
    group_input.add_argument(
        "-p",
        "--previous",
        help="File listing the locations of previously-downloaded files to \
            copy from, instead of re-downloading.",
        type=Path
    )
    group_input.add_argument(
        "-r",
        "--use-restricted",
        help="Include genomes marked as restricted (default: don't include)",
        default=False,
        action="store_true",
        dest="use_restricted"
    )
    group_input.add_argument(
        "-e",
        "--exclude-list",
        help="A text file with project names to exclude (one project per line).",
        type=Path,
        dest="exclude_list"
    )
    group_input.add_argument(
        "-j",
        "--jgi-credentials",
        help="Optional file containing JGI username (first line) and password (second line). \
            If not provided, you will be prompted to enter credentials.",
        type=Path,
        dest="credentials_file"
    )

    group_output = parser.add_argument_group("output")
    group_output.add_argument(
        "-o",
        "--outputfolder",
        help="Base directory to download all files.",
        type=Path,
        default=(Path(__file__).parent / "output"),
    )
    group_output.add_argument(
        "-s",
        "--simulate",
        help="Don't download anything, just create the folder structure and \
            taxonomy file.",
        default=False,
        action="store_true"
    )

    return parser.parse_args()


# ============================================================================
# AUTHENTICATION
# ============================================================================

def JGI_login(cookie_path, credentials_file=None):
    """
    Get (physical) cookie. Needed to download files from JGI.
    
    Args:
        cookie_path: path where the cookie file should be saved
        credentials_file: optional path to a file containing username and password
    """
    if credentials_file:
        # Read credentials from file
        try:
            with open(credentials_file, 'r') as f:
                lines = f.readlines()
                if len(lines) < 2:
                    exit("Credentials file must contain at least 2 lines (username and password)")
                user = lines[0].strip()
                password = lines[1].strip()
                if user == "" or password == "":
                    exit("Username and password cannot be empty in credentials file")
        except FileNotFoundError:
            exit(f"Credentials file not found: {credentials_file}")
        except Exception as e:
            exit(f"Error reading credentials file: {e}")
    else:
        # Prompt user for credentials
        print("Please enter your JGI credentials:")
        user = input("Username: ")
        password = getpass()
        if user == "" or password == "":
            exit("Missing info for JGI login, exiting...")

    # This should be equivalent to:
    # curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=USER_NAME' --data-urlencode 'password=USER_PASSWORD' -c cookies > /dev/null
    # See: https://genome.jgi.doe.gov/portal/help/download.jsf#api
    command = []
    command.append("curl")
    command.append("--silent")
    command.append("https://signon.jgi.doe.gov/signon/create")
    command.append("--data-urlencode")
    command.append(f"login={user}")
    command.append("--data-urlencode")
    command.append(f"password={password}")
    command.append("-c")
    command.append(f"{cookie_path}")

    print(" Logging in:")
    print("--------------------------------------------")
    proc = sp.run(command, shell=False, stdout=DEVNULL)
    try:
        proc.check_returncode()
    except sp.CalledProcessError:
        print(f"Error logging in: {proc.stderr}")
        return False
    else:
        print("--------------------------------------------\n")
        return True


# ============================================================================
# DATA RETRIEVAL & PARSING
# ============================================================================

def get_JGI_genome_list(outputfolder):
    cmd = []
    cmd.append("curl")
    cmd.append("https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/download-group?flt=&seq=all&pub=all&grp=fungi&srt=released&ord=asc")
    cmd.append("-o")
    cmd.append(str(outputfolder / "MycoCosm_Genome_list.csv"))

    print(" Downloading genome list")
    print("--------------------------------------------")
    proc = sp.run(cmd, stderr=STDOUT, encoding="utf-8")
    try:
        proc.check_returncode()
    except sp.CalledProcessError:
        print(f"Error: {proc.stderr}")
        return False
    else:
        print("--------------------------------------------\n")
        return True


def read_previous_locations(previous_location_file) -> dict:
    """Reads location of files from previous runs"""
    location_previous = {}

    with open(previous_location_file) as f:
        for line in f:
            filename, location = line.strip().split("\t")
            location_previous[filename] = location

    return location_previous


def read_excluded_projects(exclude_list_file:Path) -> set:
    """
    Reads a user-provided list of projects to exclude from the download.
    One project label per line.
    """
    exclude_projects = set()

    with open(exclude_list_file) as f:
        for line in f:
            stripped = line.strip()
            if line[0] == "#" or stripped == "":
                continue
            exclude_projects.add(stripped.split()[0])

    return exclude_projects


def remove_version(label):
    """
    Try to remove assembly version from label.
    """
    suffix = label.strip().split(" ")[-1]
    if (
        (suffix[0] == "v" or suffix[0] == "V")
        and suffix[1].isdigit()
        and suffix[-1].isdigit()
    ):
        # Assume we have one of
        # v3, v1.0, v2.0, v4.0, v1.2, v2.2, v1.1, v3.0, v2, etc.
        return " ".join(label.strip().split(" ")[:-1])
    else:
        return label


def get_final_output_folder(o, lineage_set):
    """
    Assembles the path that will contain the files for a particular organism
    based on its lineage.
    """
    # Build folder structure
    output_folder = o
    if "dikarya" in lineage_set:
        output_folder = output_folder / "DIKARYA"

        if "ascomycota" in lineage_set:
            output_folder = output_folder / "ASCOMYCOTA"

            if "pezizomycotina" in lineage_set:
                output_folder = output_folder / "PEZIZOMYCOTINA"

        elif "basidiomycota" in lineage_set:
            output_folder = output_folder / "BASIDIOMYCOTA"

            if "agaricomycotina" in lineage_set:
                output_folder = output_folder / "AGARICOMYCOTINA"

    else:
        output_folder = output_folder / "no_rank"

        if "mucoromycota" in lineage_set:
            output_folder = output_folder / "MUCOROMYCOTA"

        if "zoopagomycota" in lineage_set:
            output_folder = output_folder / "ZOOPAGOMYCOTA"

        if "chytridiomycota" in lineage_set:
            output_folder = output_folder / "CHYTRIDIOMYCOTA" / "no_rank"

    # Build last folder
    if len(lineage_set & JGI_TREE_BRANCHES) == 1:
        output_folder = output_folder / list(lineage_set & JGI_TREE_BRANCHES)[0].upper()
    else:
        output_folder = output_folder / "no_rank"

    return output_folder


def read_mycocosm_csv(csvpath: Path) -> dict:
    if not csvpath.is_file():
        exit(f"Couldn't find genome list file ({csvpath}).")

    orgdict = dict()
    # MycoCosm file is downloaded with ISO-8859-1 encoding. Use this to avoid UnicodeDecodeError.
    # See https://stackoverflow.com/questions/12468179/unicodedecodeerror-utf8-codec-cant-decode-byte-0x9c
    with codecs.open(csvpath, "r", encoding="utf-8", errors="replace") as csvfile:
        # csv reader is needed because the "name" column contains items that
        # contain commas, so a simple split() is not enough
        csvreader = csv.DictReader(csvfile)
        for l in csvreader:
            fungus = JGI_Project()
            taxid = l["NCBI Taxon"]
            shortname = l["portal"]
            fullname = l["name"]

            if taxid in ALT_TAXID:
                taxid = ALT_TAXID[taxid]

            fungus.portal = shortname
            fungus.TaxId = taxid
            if shortname in MANUAL_PROJECT_NAMES:
                fungus.name = MANUAL_PROJECT_NAMES[shortname]
            else:
                fungus.name = fullname
            fungus.org_name = remove_version(fullname)

            ncbi = NCBITaxa()
            try:
                ncbi_tax_lineage_ids = ncbi.get_lineage(taxid)
            except ValueError:
                ncbi_tax_lineage_ids = []
                print(
                    f"Can't find TaxId for {taxid} ({shortname}) with ete (try using --update)."
                )

            ncbi_tax_names = ncbi.get_taxid_translator(ncbi_tax_lineage_ids)
            lineage_set = set(
                [ncbi_tax_names[tid].lower() for tid in ncbi_tax_lineage_ids]
            )
            lineage_list_str = ",".join(
                [ncbi_tax_names[tid].lower() for tid in ncbi_tax_lineage_ids]
            )

            fungus.lineage_set = lineage_set
            fungus.lineage_list = lineage_list_str

            fungus.project_path = get_final_output_folder(Path("./"), lineage_set)

            if l["is restricted"] == "Y":
                fungus.is_restricted = True

            orgdict[shortname] = fungus

    # Remove projects we don't want to include
    for d in PORTALS_TO_REMOVE:
        if d in orgdict:
            del orgdict[d]

    return orgdict


# ============================================================================
# XML PROCESSING
# ============================================================================

def download_project_xml(cookie_path: Path, target_project: str): 
    """Given a cookie and a project name, download XML with file list."""
    command = []
    command.append("curl")
    command.append("--retry")
    command.append("5")
    command.append(f"https://genome-downloads.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism={target_project}")
    command.append("-b")
    command.append(str(cookie_path))

    try:
        xml_query = sp.run(command, encoding="utf-8", capture_output=True, check=True)
    except CalledProcessError as e:
        print(f"Command {e.cmd} failed")
        print(e.output)
        exit(f"Error while downloading XML list of files for {target_project}")

    return etree.fromstring(xml_query.stdout)


def get_JGI_xml(o: Path, cookie_path: Path, project_dict: dict):
    """
    Downloads an XML list of all files per project, for all projects.
    
    It will download data for all projects (including restricted).
    If already present, it will only try to download missing data.

    See instructions in:
    https://genome.jgi.doe.gov/portal/help/download.jsf#/api
    """
    if not cookie_path.is_file():
        exit("Error (get_JGI_xml): cookie not a valid file.")

    xmlfile = o / "MycoCosm_data.xml"
    
    # See if we already have something
    preexisting_data = set()
    if xmlfile.is_file():
        data = etree.parse(str(xmlfile))
        root = data.getroot()
        for child in root:
            if child.tag == "organismDownloads":
                preexisting_data.add(child.attrib["name"])
    else:
        # Otherwise create a new XML structure
        root = etree.Element("Data", attrib={"name":"Mycocosm"})

    # Query JGI for each project's (missing) data
    for project in project_dict.keys():
        if project in preexisting_data:
            continue

        project_xml = download_project_xml(cookie_path, project)
        root.append(project_xml)

        # Not ideal but better to write the file back each time in case process is interrupted
        with open(xmlfile, "bw") as f:
            f.write(etree.tostring(root, pretty_print=True))

    return True


# ============================================================================
# EXTRACT FILE PATHS FOR DOWNLOADING
# ============================================================================

def annotate_projects(organisms_csv: dict, xml_file: Path, hardcoded_gff_files: dict, outputfolder: Path):
    """
    Parses the XML file and extracts the genome assembly, and
    GFF3 file paths for downloading for each project.
    """
    # Auxiliary files
    gff_filenames = defaultdict(list) # holds all found GFF3 annotation files
    skipped_gffs = set()

    # Parse XML
    root = etree.parse(str(xml_file)).getroot()
    for portal_node in root:
        if portal_node.tag != "organismDownloads":
            continue

        portal = portal_node.attrib["name"]

        for child in portal_node:
            # The other option is "Mycocosm" which has a slightly different
            # set of files (at least for Trire2) but they seem to be on tape
            if child.attrib["name"] == "Files":
                Fungi_files = child
                break
        else:
            print(f"No 'Files' node in portal {portal}, skipping...")
            continue

        # Initialize variables to None to avoid UnboundLocalError
        asm_unmasked_node = None
        asm_masked_node = None
        gff_node = None
        
        for folder in Fungi_files:
            if folder.attrib["name"] == "Assembly":
                Assembly_folder = folder
                for subfolder in Assembly_folder:
                    if subfolder.attrib["name"] == "Genome Assembly (unmasked)":
                        asm_unmasked_node = subfolder
                    if subfolder.attrib["name"] == "Genome Assembly (masked)":
                        asm_masked_node = subfolder    
            if folder.attrib["name"] == "Annotation":
                Annotation_folder = folder
                for subfolder in Annotation_folder:
                    if subfolder.attrib["name"] == 'Filtered Models ("best")':
                        Filtered = subfolder
                        for subsubfolder in Filtered:
                            if subsubfolder.attrib["name"] == "Genes":
                                gff_node = subsubfolder

        # Validate that required nodes were found
        if asm_unmasked_node is None:
            raise ValueError(
                f"Portal '{portal}': Could not find 'Genome Assembly (unmasked)' node. "
                "The portal structure may have changed or this portal lacks the required assembly."
            )
        
        annotate_assembly(organisms_csv, asm_unmasked_node)

        # A few projects don't have unmasked assemblies. Assign masked ones in this case.
        if not organisms_csv[portal].assembly_file:
            print(f"Portal {portal} missing unmasked assembly. ", end="")
            if asm_masked_node is None:
                print("Could not find masked assembly as fallback. Skipping assembly annotation for this portal.")
            else:
                print("Attempting masked version...")
                annotate_missing(organisms_csv, asm_masked_node)
        
        # Validate GFF node before annotation
        if gff_node is None:
            print(f"Portal '{portal}': Warning - Could not find GFF node (Filtered Models > best > Genes). "
                  "Skipping GFF annotation for this portal.")
        else:
            annotate_gff(organisms_csv, gff_node, hardcoded_gff_files, gff_filenames, skipped_gffs)
        
    # Write down all found annotation files.
    repeated_file = outputfolder / "List_gene_gff_filenames.txt"
    with codecs.open(repeated_file, "w", encoding="utf-8") as f:
        for portal in sorted(gff_filenames):
            f.write(f"{portal} ({organisms_csv[portal].name})\n")
            f.write(f"\t{organisms_csv[portal].project_path}\n")
            for filename, dt in gff_filenames[portal]:
                if filename in skipped_gffs:
                    f.write(f"\t{dt.strftime('%Y-%m-%d')}\t{filename} (SKIPPED)")
                else:
                    f.write(f"\t{dt.strftime('%Y-%m-%d')}\t{filename}")

                if filename == organisms_csv[portal].gff_file:
                    f.write(" *\n")
                else:
                    f.write("\n")
            f.write("\n")
        
    # Report if we're missing annotation files
    portals_missing_gff_file = [
        p.portal for p in organisms_csv.values() if p.gff_file == ""
    ]
    if len(portals_missing_gff_file) > 0:
        missing_gff_string = ", ".join(portals_missing_gff_file)
        print(f"Missing gffs from: {missing_gff_string}")


def annotate_assembly(organisms_csv, xml_unmasked_assembly):
    duplicated_names = []
    portal = None  # Initialize to None to avoid UnboundLocalError if iterator yields no elements
    # Projects' names differ between what's annotated in the CSV and XML files.
    #
    # For example:
    # CSV: 'Mortierella humilis PMI_1414 v1.0'
    # XML: 'Mortierella humilis PMI_1414 - Glomeribacter phylotype 2 Fungal Standard Draft'
    for element in xml_unmasked_assembly.iterdescendants():
        filename = element.attrib["filename"]
        url = element.attrib["url"]
        portal = url.split("/")[2]
        file_size = int(element.attrib["sizeInBytes"])

        if (
            "MitoAssembly" in filename
            or "MitoScaffolds" in filename
            or "PrimaryAssemblyScaffolds" in filename
            or "SecondaryAssemblyScaffolds" in filename
            or filename in EXCLUDE_ASSEMBLIES
            or portal in PORTALS_TO_REMOVE
            or filename.endswith("txt")
        ):
            continue

        try:
            organisms_csv[portal].assembly_file = filename
            organisms_csv[portal].assembly_url = url
            organisms_csv[portal].assembly_size = file_size
        except KeyError:
            print(
                f"WARNING! Portal {portal} from XML file not found in CSV list."
            )

        # Add filenames to list of found assembly files for this portal to check for multiple files
        duplicated_names.append(filename)

    # Check to see if there are multiple assembly files for this portal
    if len(duplicated_names) > 1:
        if portal is not None:
            print(f"WARNING! Portal {portal} has more than one assembly file.")     
            dup_names = ", ".join(duplicated_names)
            print(f"{portal}\t{dup_names}")
            print("")
        else:
            print("WARNING! Multiple assembly files found but portal information is unavailable.")
    elif len(duplicated_names) == 0:
        if portal is not None:
            print(f"Portal {portal} has no (unmasked) assembly file.")
        else:
            print("No (unmasked) assembly files found in the provided XML data.")

    return


def annotate_missing(organisms_csv, xml_masked_assembly):
    duplicated_names = []
    portal = None  # Initialize to None to avoid UnboundLocalError if iterator yields no elements
    for element in xml_masked_assembly.iterdescendants():
        filename = element.attrib["filename"]
        url = element.attrib["url"]
        portal = url.split("/")[2]
        file_size = int(element.attrib["sizeInBytes"])

        if (
            "MitoAssembly" in filename
            or "MitoScaffolds" in filename
            or "PrimaryAssemblyScaffolds" in filename
            or "SecondaryAssemblyScaffolds" in filename
            or filename in EXCLUDE_ASSEMBLIES
            or portal in PORTALS_TO_REMOVE
            or filename.endswith("txt")
        ):
            continue

        organisms_csv[portal].assembly_file = filename
        organisms_csv[portal].assembly_url = url
        organisms_csv[portal].assembly_size = file_size

        duplicated_names.append(filename)

    if len(duplicated_names) > 1:
        if portal is not None:
            print(f"Warning: portal {portal} has more than one masked assembly file!")     
            dup_names = ", ".join(duplicated_names)
            print(f"{portal}\t{dup_names}")
            print("")
        else:
            print("Warning: Multiple masked assembly files found but portal information is unavailable.")
    elif len(duplicated_names) == 0:
        if portal is not None:
            print(f"Portal {portal} has no (masked) assembly file")
        else:
            print("No (masked) assembly files found in the provided XML data.")

    return


def get_hardcoded_gffs():
    hgfs = dict()
    tsv_file = Path(__file__).parent / "hardcoded_gff_files.tsv"
    if not tsv_file.is_file():
        print("Missing 'hardcoded_gff_files.tsv'.")
        return hgfs

    with open(tsv_file) as f:
        for line in f:
            if line[0] == "#":
                continue
            portal, gff_file = line.strip().split("\t")
            hgfs[portal] = gff_file

    return hgfs


def annotate_gff(organisms_csv, gff_node, hardcoded_gff_files, gff_filenames, skipped_gffs):
    for element in gff_node.iterdescendants():
        filename = element.attrib["filename"]
        url = element.attrib["url"]
        portal = url.split("/")[2]
        timestamp = element.attrib["timestamp"]
        file_size = int(element.attrib["sizeInBytes"])

        if portal in PORTALS_TO_REMOVE:
            continue

        org = organisms_csv[portal]

        # Convert the timestamp-date into a datetime object
        # For example: 'Sun Oct 12 11:02:03 PDT 2014'
        wday, month, day, HMS, tz, y = timestamp.split(" ")
        stringtimestamp = f"{US_DAY_NUMBERS[wday]}"
        stringtimestamp += f" {US_MONTH_NUMBERS[month]}"
        stringtimestamp += f" {day}"
        stringtimestamp += f" {HMS}"
        stringtimestamp += f" {US_TIMEZONES[tz]}"
        stringtimestamp += f" {y}"
        dt_timestamp = datetime.strptime(stringtimestamp, "%w %m %d %H:%M:%S %z %Y")

        gff_filenames[portal].append((filename, dt_timestamp))

        # Conditions for skpping a gff file
        
        # Condition 1: if we have a hardcoded gff file for this portal, accept ones that have necessary data
        if portal in hardcoded_gff_files:
            if hardcoded_gff_files[portal] == filename:
                org.gff_file = filename
                org.gff_url = url
                org.gff_timestamp = dt_timestamp
                org.gff_size = file_size
                continue
            else:
                skipped_gffs.add(filename)
                continue
        
        # Condition 2: if the filename is in the list of unwanted gff files, skip
        if filename in IGNORE_GFF_FILENAMES:
            skipped_gffs.add(filename)
            continue

        # Condition 3: if the filename contains certain keywords, skip
        if "proteins" in filename.lower() \
            or "secondary_alleles" in filename.lower() \
            or "promoter_regions" in filename.lower():
            skipped_gffs.add(filename)
            continue

        # Condition 4: if the filename doesn't end with "gff.gz" or "gff3.gz", skip
        if filename[-6:] == "gtf.gz" or filename[-3:] == "tgz" or filename[-2:] != "gz":
            skipped_gffs.add(filename)
            continue

        if org.gff_file == "":
            # First filename for this portal, record
            org.gff_file = filename
            org.gff_url = url
            org.gff_timestamp = dt_timestamp
            org.gff_size = file_size
        else:
            if filename[-7:] == "gff3.gz" and org.gff_file[-6:] == "gff.gz":
                # Had gff, now got gff3. Update
                org.gff_file = filename
                org.gff_url = url
                org.gff_timestamp = dt_timestamp
                org.gff_size = file_size
            elif (
                filename[-7:] == "gff3.gz"
                and org.gff_file[-7:] == "gff3.gz"
                or filename[-6:] == "gff.gz"
                and org.gff_file[-6:] == "gff.gz"
            ):
                if dt_timestamp > org.gff_timestamp:
                    # Keep the newest
                    org.gff_file = filename
                    org.gff_url = url
                    org.gff_timestamp = dt_timestamp
                    org.gff_size = file_size
            elif filename[-6:] == "gff.gz" and org.gff_file[-7:] == "gff3.gz":
                # A gff file appeared but we had gff3, skip
                continue
            else:
                print(f" annotate_gff unexpected case ({portal})")
                print(org.gff_file)
                print(filename)
                print()
                continue

    return


# ============================================================================
# FILE OPERATIONS
# ============================================================================

def download_file(url, local_path, cookie):
    """
    Handles downloads for files
    """
    command = [
        "curl",
        f"https://genome.jgi.doe.gov{url}",
        "-b",
        f"{cookie}",
        "--retry",
        "5",
        "--output",
        f"{local_path}",
    ]
    try:
        sp.run(command, shell=False, stdout=DEVNULL, stderr=STDOUT)
    except CalledProcessError as e:
        print(f" Error: {str(e.output)}")
        return False
    else:
        return True


def get_aux_files(args, cookie_path) -> None:
    if args.update:
        print("Updating NCBI taxonomy database")
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()
        print("  Done")

    if args.getgenomelist or args.getxml:
        #  Output folder
        o = args.outputfolder
        if not o.is_dir():
            o.mkdir(exist_ok=True, parents=True)

    if args.getgenomelist:
        print("Downloading genome list (.csv)")
        if get_JGI_genome_list(o):
            print("  Done")
        else:
            print("  Cannot download MycoCosm's genome list...")

    if args.getxml:
        mycocosm_csv = o / "MycoCosm_Genome_list.csv"
        print("Downloading XML file")
        if not JGI_login(cookie_path, args.credentials_file):
            exit("Error: Cannot log in to JGI...")
        elif not mycocosm_csv.is_file():
            exit("Error downloading XML files: download genome list first")

        if get_JGI_xml(o, cookie_path, read_mycocosm_csv(mycocosm_csv)):
            print("  Done")
        else:
            print("  Cannot download MycoCosm's XML file...")

    if args.getprevious:
        print("Downloading location list of previously downloaded files")
        p = args.getprevious
        if not p.is_dir():
            exit("Error: given base folder is not a valid folder")
        with open(
            Path(__file__).parent / "previously_downloaded_files.tsv", "w"
        ) as f:
            for gz in p.glob("**/*.gz"):
                f.write(f"{gz.name}\t{gz.parent}\n")
        print("  Done")
    return


def parse_xml(xml):
    """
    Reads the XML file and returns the nodes containing portals' main data
    """
    doc = etree.parse(str(xml))
    root = doc.getroot()
    for child in root:
        if child.attrib["name"] == "Files":
            Fungi_files = child
            break

    for folder in Fungi_files:
        if folder.attrib["name"] == "Assembly":
            for subfolder in folder:
                if subfolder.attrib["name"] == "Genome Assembly (unmasked)":
                    asm_unmasked_node = subfolder
                if subfolder.attrib["name"] == "Genome Assembly (masked)":
                    asm_masked_node = subfolder
                    

        if folder.attrib["name"] == "Annotation":
            Annotation_folder = folder
            for subfolder in Annotation_folder:
                if subfolder.attrib["name"] == 'Filtered Models ("best")':
                    Filtered = subfolder
                    for subsubfolder in Filtered:
                        if subsubfolder.attrib["name"] == "Genes":
                            gff_node = subsubfolder

    return asm_unmasked_node, asm_masked_node, gff_node


# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def main():
    # ========== INITIALIZATION ==========
    args = command_parser()
    cookie_path = Path(__file__).parent / "jgi_cookie.txt"

    # ========== AUXILIARY OPERATIONS ==========
    # Handle --update, --getgenomelist, --getxml, --getprevious options and exit
    if args.update or args.getgenomelist or args.getxml or args.getprevious:
        get_aux_files(args, cookie_path)
        exit("All data fetched")
    
    # ========== INPUT VALIDATION ==========
    if args.xml is None or not args.xml.is_file():
        exit("Error: invalid file specified with --xml")
    
    if args.csv is None or not args.csv.is_file():
        exit("Error: invalid file specified with --csv")

    # ========== SETUP & DATA LOADING ==========
    # Create output directory
    o = args.outputfolder
    if not o.is_dir():
        o.mkdir(exist_ok=True, parents=True)  # recursive folder creation

    # Load genome metadata from CSV
    organisms_csv = dict()
    organisms_csv = read_mycocosm_csv(args.csv)
    print(f"Total number of genomes: {len(organisms_csv)}\n")

    # Load optional filter lists
    location_previous = dict()
    if args.previous:
        location_previous = read_previous_locations(args.previous)

    exclude_projects = {}
    if args.exclude_list:
        exclude_projects = read_excluded_projects(args.exclude_list)

    # ========== PARSE & ANNOTATE METADATA ==========
    # Populate all the info in each JGI_Project object from XML
    hardcoded_gff_files = get_hardcoded_gffs()
    annotate_projects(organisms_csv, args.xml, hardcoded_gff_files, o)

    # ========== LOGIN & DOWNLOAD PREPARATION ==========
    if args.simulate:
        print("\nBeginning simulation")
    else:
        print("\nBeginning file download")
        if not JGI_login(cookie_path, args.credentials_file):
            exit("Error: Failed to log in to JGI...")

    # ========== GENERATE DOWNLOAD LIST & PROCESS FILES ==========
    # Create a checkpoint file listing all files for each portal
    download_list_file = f"JGI_download_list_{time.strftime('%Y-%m-%d', time.localtime())}.txt"
    with open(o / download_list_file, "w", encoding="utf-8") as f:
        for portal in organisms_csv:
            fungus = organisms_csv[portal]
            f.write(f"{portal} ({fungus.name})\n")
            f.write(f"Assembly:\t{fungus.assembly_file}\n")
            f.write(f"GFF:\t\t{fungus.gff_file}\n\n")

    # Finally, get the files
    needed = 0
    copied = 0
    downloaded = 0
    pre_existing = 0
    in_previous = 0
    with open(o / "JGI_taxonomy.tsv", "w", encoding="utf-8") as tf:
        tf.write(
            "Short name\tAccession\tTaxId\tName\tPath\tAssembly file\tGFF file\tlineage\n"
        )

        for portal, fungus in organisms_csv.items():
            if not args.use_restricted and fungus.is_restricted:
                continue
            if portal in exclude_projects:
                continue
            
            if fungus.assembly_file == "" or fungus.gff_file == "":
                if fungus.assembly_file == "" and fungus.gff_file == "":
                    missing_string = "(assembly and gff)"
                elif fungus.assembly_file == "":
                    missing_string = "(assembly)"
                else:
                    missing_string = "(gff)"
                print(f" Warning! Missing file for {portal}. ", end="")
                print(f" {missing_string}. See {download_list_file}")
                continue

            needed += 2

            base_folder = Path(fungus.project_path) / portal
            output_folder = o / base_folder
            if not output_folder.is_dir():
                output_folder.mkdir(exist_ok=True, parents=True)

            tf.write(portal)
            tf.write(f"\t{portal}")
            tf.write(f"\t{fungus.TaxId}")
            tf.write(f"\t{fungus.name}")
            tf.write(f"\t{base_folder}")
            tf.write(f"\t{fungus.assembly_file}")
            tf.write(f"\t{fungus.gff_file}")
            tf.write(f"\t{fungus.lineage_list}\n")

            asm = output_folder / fungus.assembly_file
            gff = output_folder / fungus.gff_file

            # Assembly
            if asm.is_file() and asm.stat().st_size > 0.9 * fungus.assembly_size:
                # File already there, skipping
                pre_existing += 1
            elif fungus.assembly_file in location_previous:
                if args.simulate:
                    in_previous += 1
                else:
                    old_file = (
                        Path(location_previous[fungus.assembly_file])
                        / fungus.assembly_file
                    )
                    try:
                        shutil.copy(old_file, output_folder)
                    except:
                        exit(f"Error: cannot copy {old_file} into {output_folder}")
                    else:
                        copied += 1
            else:
                if not args.simulate:
                    if download_file(fungus.assembly_url, asm, cookie_path):
                        downloaded += 1
                    else:
                        print(f"Warning: could not download assembly {asm}")

            # GFF
            if gff.is_file() and gff.stat().st_size > 0.9 * fungus.gff_size:
                pre_existing += 1
            elif fungus.gff_file in location_previous:
                if args.simulate:
                    in_previous += 1
                else:
                    old_file = (
                        Path(location_previous[fungus.gff_file]) / fungus.gff_file
                    )
                    try:
                        shutil.copy(old_file, output_folder)
                    except:
                        exit(f"Error: cannot copy {old_file} into {output_folder}")
                    else:
                        copied += 1
            else:
                if not args.simulate:
                    if download_file(fungus.gff_url, gff, cookie_path):
                        downloaded += 1
                    else:
                        print(f"Warning: could not download gff {gff}")

    print("...done\n")

    # ========== SUMMARY ==========
    print(f"Portals: {len(organisms_csv)}")
    print(f"Files needed: {needed}")
    total_got = copied + downloaded + pre_existing + in_previous
    print(f"Got {copied} (copied)", end="")
    print(f" + {downloaded} (downloaded)", end="")
    print(f" + {pre_existing} (pre-existing)", end="")
    print(f" + {in_previous} (in-previous)", end="")
    print(f" = {total_got}")


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    main()
