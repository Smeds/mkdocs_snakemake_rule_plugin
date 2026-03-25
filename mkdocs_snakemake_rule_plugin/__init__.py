
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("mkdocs-snakemake-rule-plugin")
except PackageNotFoundError:
    __version__ = "unknown"
