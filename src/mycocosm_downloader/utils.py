def remove_version(label: str) -> str:
    """Remove assembly version from label if present."""
    parts = label.strip().split(" ")
    suffix = parts[-1]
    if (suffix[0].lower() == "v" and suffix[1].isdigit() and suffix[-1].isdigit()):
        return " ".join(parts[:-1])
    return label