# MkDocs Snakemake Rule

MkDocs Plugin used to extract rule information.

To use this plugin, install it with pip in the same environment as MkDocs:

```
pip install mkdocs-snakemake-rule-plugin
```

Then add the following entry to the MkDocs config file:

```yml
plugins:
  - snakemake-rule:
      rule_folders:
        - 'workflow/rules'
```

In your target file, add a tag to be replaced with format

SNAKEMAKE_RULE_SOURCE__filename__rulename

- filename is where the rule is stored, can be without '.skm'
- rulename is the rule that will be extract


```
#SNAKEMAKE_RULE_SOURCE__fastp__fastp_pe#
```
