from dataclasses import dataclass, field
from pathlib import Path
from typing import Set, List


@dataclass
class MycoCosmProject:
    """Store information about each MycoCosm sequencing project (portal)."""

    taxid: str = ""
    lineage_set: Set[str] = field(default_factory=set)
    lineage_list: List[str] = field(default_factory=list)

    portal: str = ""
    name: str = ""
    org_name: str = ""

    assembly_file: str = ""
    assembly_url: str = ""
    assembly_size: int = 0

    gff_file: str = ""
    gff_url: str = ""
    gff_timestamp: str = ""
    gff_size: int = 0

    protein_file: str = ""
    protein_url: str = ""
    protein_size: int = 0

    project_path: Path | None = None

    date: str = ""   # later upgrade to datetime
    is_restricted: bool = False
