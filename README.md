# MkDocs Snakemake Rule

MkDocs Plugin used to extract rule information.

## Installation

To use this plugin, install it with pip in the same environment as MkDocs:

```
pip install mkdocs-snakemake-rule-plugin
```

## Config

Then add the following entry to the MkDocs config file:

```yml
plugins:
  - snakemake-rule:
      rule_folders:
        - 'workflow/rules'
      schemas:
        - 'workflow/schemas/rules.schema.yaml'
```

## Tag 

### format

Format of tags are:
- `#SNAKEMAKE_RULE_SOURCE__filename__rulename#` : for source code extraction
- `#  SNAKEMAKE_RULE_TABLE__filename__rulename#` : for table generation

Tab parts:
- first section identifies if rule source or table should be created
- filename is where the rule is stored, can be without '.skm'
- rulename is the rule that will be extract


### Usage

Add tags to your target file and they will be replaced with rule source or tables

MkDocs Markdown example
```
# Rule information
## Source 
#SNAKEMAKE_RULE_SOURCE__fastp__fastp_pe#

## Parameters
#SNAKEMAKE_RULE_TABLEE__fastp__fastp_pe#

```
