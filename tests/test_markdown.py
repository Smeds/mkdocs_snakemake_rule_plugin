import warnings

import pytest
import yaml

from mkdocs_snakemake_rule_plugin.markdown import (
    IncorrectTagError,
    RuleMissingSchemaFileError,
    extract_snakemake_rule,
    extract_snakemake_rule_section,
    markdown_gen,
    markdown_source,
    markdown_table,
    parse_variable,
    remove_brackets,
    remove_comment,
    remove_indent,
    remove_temp_and_output,
    remove_trailing_empty_lines,
    replace_newline,
)

# ── Fixtures ──────────────────────────────────────────────────────────────────

SAMPLE_SMK = """\
rule align_reads:
    input:
        reads = "data/{sample}.fastq",
        ref = "reference/genome.fa"
    output:
        bam = "results/{sample}.bam"
    params:
        extra = "--very-sensitive"
    log:
        "logs/{sample}.log"
    shell:
        "bwa mem {params.extra} {input.ref} {input.reads} > {output.bam}"

checkpoint split_file:
    input:
        infile = "data/input.txt"
    output:
        outdir = directory("results/split/")
    shell:
        "split.sh {input} {output}"

rule second_rule:
    input:
        "data/file.txt"
    output:
        "results/file.txt"
    shell:
        "cp {input} {output}"

rule rule_with_unpack:
    input:
        unpack(lambda wildcards:  get_input_haplotagged_bam(wildcards, config, set_type="T")),
        vntr=config.get("severus_t_only", {}).get("vntr", ""),
        pon=config.get("severus_t_only", {}).get("pon", ""),
"""

SAMPLE_SCHEMA = {
    "properties": {
        "align_reads": {
            "properties": {
                "input": {
                    "properties": {
                        "reads": {"description": "Input reads FASTQ file"},
                        "ref": {"description": "Reference genome FASTA file"},
                    }
                },
                "output": {
                    "properties": {
                        "bam": {"description": "Output BAM file"},
                    }
                },
            }
        }
    }
}


@pytest.fixture
def smk_file(tmp_path):
    path = tmp_path / "sample.smk"
    path.write_text(SAMPLE_SMK)
    return str(path)


@pytest.fixture
def schema_file(tmp_path):
    path = tmp_path / "schema.yaml"
    path.write_text(yaml.dump(SAMPLE_SCHEMA))
    return str(path)


@pytest.fixture
def rule_folder(tmp_path):
    rules_dir = tmp_path / "rules"
    rules_dir.mkdir()
    (rules_dir / "sample.smk").write_text(SAMPLE_SMK)
    return str(rules_dir)


@pytest.fixture
def configured_gen(rule_folder, schema_file):
    gen = markdown_gen()
    gen.set_config({"rule_folders": [rule_folder], "schemas": [schema_file]})
    return gen


# ── Pure utility functions ─────────────────────────────────────────────────────

class TestReplaceNewline:
    def test_replaces_single_newline(self):
        assert replace_newline("hello\nworld") == "hello<br />world"

    def test_no_newline_unchanged(self):
        assert replace_newline("hello world") == "hello world"

    def test_multiple_newlines(self):
        assert replace_newline("a\nb\nc") == "a<br />b<br />c"


class TestRemoveComment:
    def test_no_op_on_normal_string(self):
        # remove_comment does a literal string replace of r"#.*", not a regex substitution
        assert remove_comment("no comment here") == "no comment here"

    def test_removes_literal_pattern(self):
        # Only the literal string "#.*" is removed, not real inline comments
        assert remove_comment("value #.*") == "value "

    def test_real_inline_comment_not_removed(self):
        # A real inline comment is NOT stripped (known limitation of the implementation)
        assert remove_comment("value # comment") == "value # comment"


class TestRemoveIndent:
    def test_multiple_spaces_collapsed(self):
        assert remove_indent("a  b") == "a b"

    def test_multiple_tabs_collapsed(self):
        assert remove_indent("a\t\tb") == "a b"

    def test_single_space_unchanged(self):
        assert remove_indent("a b") == "a b"

    def test_leading_spaces_collapsed(self):
        assert remove_indent("   value") == " value"


class TestRemoveBrackets:
    def test_removes_leading_and_trailing_brackets(self):
        assert remove_brackets("[ value1, value2 ]") == " value1, value2 "

    def test_no_brackets_unchanged(self):
        assert remove_brackets("value") == "value"

    def test_string_starting_with_bracket_stripped(self):
        assert remove_brackets("[value]") == "value"


class TestRemoveTempAndOutput:
    def test_removes_temp_wrapper(self):
        assert remove_temp_and_output("temp(results/{sample}.bam)") == "results/{sample}.bam"

    def test_removes_directory_wrapper(self):
        assert remove_temp_and_output("directory(results/split/)") == "results/split/"

    def test_no_wrapper_unchanged(self):
        assert remove_temp_and_output("results/{sample}.bam") == "results/{sample}.bam"

    def test_temp_with_leading_spaces(self):
        assert remove_temp_and_output("   temp(results/{sample}.bam)") == "results/{sample}.bam"


class TestParseVariable:
    def test_balanced_parentheses_returns_immediately(self):
        rows = iter([])
        assert parse_variable("expand('{sample}')", rows) == "expand('{sample}')"

    def test_unbalanced_reads_next_line(self):
        rows = iter(["second_part)"])
        result = parse_variable("expand('{sample}',", rows)
        assert result == "expand('{sample}',second_part)"

    def test_balanced_brackets(self):
        rows = iter([])
        assert parse_variable("[a, b]", rows) == "[a, b]"

    def test_nested_unbalanced_reads_multiple_lines(self):
        rows = iter(["middle,", "end))"])
        result = parse_variable("expand(('{s}',", rows)
        assert result == "expand(('{s}',middle,end))"


class TestMarkdownSource:
    def test_wraps_in_fenced_code_block(self):
        result = markdown_source("rule foo:\n    pass")
        assert result == "```\nrule foo:\n    pass\n```"

    def test_empty_string(self):
        assert markdown_source("") == "```\n\n```"


class TestRemoveTrailingEmptyLines:
    def test_strips_trailing_newline(self):
        assert remove_trailing_empty_lines("rule foo:\n    pass\n") == "rule foo:\n    pass"

    def test_strips_multiple_trailing_newlines(self):
        assert remove_trailing_empty_lines("rule foo:\n    pass\n\n\n") == "rule foo:\n    pass"

    def test_no_trailing_newline_unchanged(self):
        assert remove_trailing_empty_lines("rule foo:\n    pass") == "rule foo:\n    pass"


class TestExtractSnakemakeRuleSection:
    def test_empty_parts_returns_data(self):
        assert extract_snakemake_rule_section([], {"key": "value"}) == {"key": "value"}

    def test_non_empty_parts_still_returns_data(self):
        # The function recurses until parts is empty but ignores part values
        assert extract_snakemake_rule_section(["a", "b"], {"key": "value"}) == {"key": "value"}


# ── extract_snakemake_rule ─────────────────────────────────────────────────────

class TestExtractSnakemakeRule:
    def test_extracts_named_rule(self, smk_file):
        result = extract_snakemake_rule(smk_file, "align_reads")
        assert "rule align_reads:" in result
        assert 'reads = "data/{sample}.fastq"' in result
        assert 'bam = "results/{sample}.bam"' in result

    def test_extracts_checkpoint(self, smk_file):
        result = extract_snakemake_rule(smk_file, "split_file")
        assert "checkpoint split_file:" in result
        assert 'infile = "data/input.txt"' in result

    def test_missing_rule_returns_empty_string(self, smk_file):
        result = extract_snakemake_rule(smk_file, "nonexistent_rule")
        assert result == ""

    def test_rule_content_stops_before_next_rule(self, smk_file):
        result = extract_snakemake_rule(smk_file, "align_reads")
        assert "second_rule" not in result
        assert "split_file" not in result

    def test_extracts_last_rule_in_file(self, smk_file):
        result = extract_snakemake_rule(smk_file, "second_rule")
        assert "rule second_rule:" in result
        assert 'cp {input} {output}' in result

    def test_extracts_rule_with_unpack(self, smk_file):
        result = extract_snakemake_rule(smk_file, "rule_with_unpack")
        assert "rule rule_with_unpack:" in result


# ── markdown_table ─────────────────────────────────────────────────────────────

class TestMarkdownTable:
    RULE_SOURCE = (
        "rule align_reads:\n"
        "    input:\n"
        '        reads = "data/{sample}.fastq",\n'
        '        ref = "reference/genome.fa"\n'
        "    output:\n"
        '        bam = "results/{sample}.bam"\n'
    )

    def test_generates_table_header(self):
        schema = {
            "properties": {
                "input": {"properties": {"reads": {"description": "Input reads"}}}
            }
        }
        result = markdown_table(self.RULE_SOURCE, schema)
        assert "| Rule parameters | Key | Value | Description |" in result
        assert "| --- | --- | --- | --- |" in result

    def test_includes_variable_key_and_description(self):
        schema = {
            "properties": {
                "input": {"properties": {"reads": {"description": "Input reads FASTQ file"}}}
            }
        }
        result = markdown_table(self.RULE_SOURCE, schema)
        assert "reads" in result
        assert "Input reads FASTQ file" in result

    def test_multiple_sections(self):
        schema = {
            "properties": {
                "input": {"properties": {"reads": {"description": "Reads"}}},
                "output": {"properties": {"bam": {"description": "Output BAM"}}},
            }
        }
        result = markdown_table(self.RULE_SOURCE, schema)
        assert "input" in result
        assert "output" in result
        assert "Output BAM" in result

    def test_schema_value_overrides_extracted_value(self):
        schema = {
            "properties": {
                "input": {
                    "properties": {
                        "reads": {"description": "Reads", "value": "fixed_value"}
                    }
                }
            }
        }
        result = markdown_table(self.RULE_SOURCE, schema)
        assert "fixed_value" in result

    def test_rule_with_unpack(self, smk_file):
        rule = extract_snakemake_rule(smk_file, "rule_with_unpack")
        schema = {
            "properties": {
                "input": {
                    "properties": {
                        "unpack": {"description": "Reads"},
                        "vntr": {"description": "vtr"},
                        "pon": {"description": "pon"}
                    }
                }
            }
        }
        result = markdown_table(rule, schema)
        print(result)
        assert '`lambda wildcards: get_input_haplotagged_bam(wildcards, config, set_type="T")`' in result


class TestMarkdownGenSetConfig:
    def test_warns_when_no_rule_folders(self):
        gen = markdown_gen()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            gen.set_config({})
        messages = [str(warning.message) for warning in w]
        assert any("rule folders" in m.lower() for m in messages)

    def test_warns_when_no_schemas(self):
        gen = markdown_gen()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            gen.set_config({})
        messages = [str(warning.message) for warning in w]
        assert any("schema" in m.lower() for m in messages)

    def test_sets_rule_folders(self, rule_folder):
        gen = markdown_gen()
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            gen.set_config({"rule_folders": [rule_folder]})
        assert rule_folder in gen.config_rule_folders

    def test_loads_schema_file(self, rule_folder, schema_file):
        gen = markdown_gen()
        gen.set_config({"rule_folders": [rule_folder], "schemas": [schema_file]})
        assert len(gen.config_schemas) == 1
        assert "properties" in gen.config_schemas[0]

    def test_initialises_empty_extracted_rules(self, rule_folder, schema_file):
        gen = markdown_gen()
        gen.set_config({"rule_folders": [rule_folder], "schemas": [schema_file]})
        assert gen.config_extracted_rules == {}


class TestMarkdownGenFindFile:
    def test_finds_file_with_explicit_smk_extension(self, configured_gen):
        paths = configured_gen.find_file("sample.smk")
        assert len(paths) == 1
        assert paths[0].endswith("sample.smk")

    def test_finds_file_without_extension(self, configured_gen):
        paths = configured_gen.find_file("sample")
        assert len(paths) == 1
        assert paths[0].endswith("sample.smk")

    def test_returns_empty_list_for_missing_file(self, configured_gen):
        paths = configured_gen.find_file("nonexistent")
        assert paths == []

    def test_none_returns_all_smk_files(self, configured_gen):
        paths = configured_gen.find_file(None)
        assert len(paths) >= 1
        assert all(p.endswith(".smk") for p in paths)


# ── markdown_gen.extract_rule_source ──────────────────────────────────────────

class TestMarkdownGenExtractRuleSource:
    def test_extracts_rule_by_name(self, configured_gen):
        result = configured_gen.extract_rule_source("sample.smk", "align_reads")
        assert result is not None
        assert "rule align_reads:" in result["source"]
        assert "file_path" in result

    def test_caches_extracted_rule(self, configured_gen):
        configured_gen.extract_rule_source("sample.smk", "align_reads")
        assert "align_reads" in configured_gen.config_extracted_rules

    def test_second_call_returns_cached_object(self, configured_gen):
        result1 = configured_gen.extract_rule_source("sample.smk", "align_reads")
        result2 = configured_gen.extract_rule_source("sample.smk", "align_reads")
        assert result1 is result2

    def test_file_path_recorded_in_cache(self, configured_gen, rule_folder):
        result = configured_gen.extract_rule_source("sample.smk", "align_reads")
        assert rule_folder in result["file_path"]


# ── markdown_gen.get_markdown ─────────────────────────────────────────────────

class TestMarkdownGenGetMarkdown:
    def test_markdown_without_tags_unchanged(self, configured_gen):
        md = "# Hello\n\nSome text with no tags."
        assert configured_gen.get_markdown(md) == md

    def test_source_tag_replaced_with_code_block(self, configured_gen):
        md = "#SNAKEMAKE_RULE_SOURCE__sample__align_reads#"
        result = configured_gen.get_markdown(md)
        assert "```" in result
        assert "rule align_reads:" in result
        assert "#SNAKEMAKE_RULE_SOURCE__sample__align_reads#" not in result

    def test_table_tag_replaced_with_table(self, configured_gen):
        md = "#SNAKEMAKE_RULE_TABLE__sample__align_reads#"
        result = configured_gen.get_markdown(md)
        assert "| Rule parameters |" in result
        assert "#SNAKEMAKE_RULE_TABLE__sample__align_reads#" not in result

    def test_tag_embedded_in_surrounding_text(self, configured_gen):
        md = "See rule: #SNAKEMAKE_RULE_SOURCE__sample__align_reads# for details."
        result = configured_gen.get_markdown(md)
        assert "See rule:" in result
        assert "for details." in result
        assert "rule align_reads:" in result

    def test_incorrect_source_tag_raises(self, configured_gen):
        md = "#SNAKEMAKE_RULE_SOURCE__onlytwoparts#"
        with pytest.raises(IncorrectTagError):
            configured_gen.get_markdown(md)

    def test_incorrect_table_tag_raises(self, configured_gen):
        md = "#SNAKEMAKE_RULE_TABLE__onlytwoparts#"
        with pytest.raises(IncorrectTagError):
            configured_gen.get_markdown(md)

    def test_table_tag_missing_schema_raises(self, rule_folder, tmp_path):
        # Schema exists but does not contain the requested rule
        schema_without_rule = {"properties": {}}
        schema_path = tmp_path / "empty_schema.yaml"
        schema_path.write_text(yaml.dump(schema_without_rule))

        gen = markdown_gen()
        gen.set_config({"rule_folders": [rule_folder], "schemas": [str(schema_path)]})

        md = "#SNAKEMAKE_RULE_TABLE__sample__align_reads#"
        with pytest.raises(RuleMissingSchemaFileError):
            gen.get_markdown(md)
