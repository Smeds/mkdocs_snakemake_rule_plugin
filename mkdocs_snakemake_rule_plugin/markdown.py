import glob
import os
import re
import yaml
import warnings

TAG_SNAKEMAKE_SOURCE_RE = r"#SNAKEMAKE_RULE_SOURCE[A-Za-z0-9_-]+#"
TAG_SNAKEMAKE_TABLE_RE = r"#SNAKEMAKE_RULE_TABLE[A-Za-z0-9_-]+#"


class IncorrectTagError(ValueError):
    'Incorrect tag'
    pass


class MissingConfigSettingError(ValueError):
    'Missing config setting'
    pass


class RuleFileMissingError(IOError):
    'Missing config setting'
    pass


class RuleMissingInSmkFileError(ValueError):
    'Missing config setting'
    pass


class RuleMissingSchemaFileError(ValueError):
    'Missing config setting'
    pass


def extract_snakemake_rule_section(parts, data):
    if len(parts) > 0:
        return extract_snakemake_rule_section(parts[1:], data)
    else:
        return data


def extract_snakemake_rule(file_path, rule):
    rule_content = ""
    with open(file_path, 'r') as reader:
        for line in reader:
            if line.startswith('rule') and rule in line:
                rule_content += line
                for line in reader:
                    if not line.startswith("rule") or line.startswith("def") or re.search(r"^[A-Za-z_-]+"):
                        rule_content += line
                    else:
                        return rule_content
    return rule_content


def markdown_source(code_item):
    return "```\n" + code_item + "\n```"


def replace_newline(value):
    return value.replace("\n", "<br />")


def remove_indent(value):
    return re.sub(r"[ ]{2,}", " ", re.sub(r"\t{2,}", " ", value))


def remove_brackets(value):
    if re.match(r"[ ]*\[", value):
        value = re.sub(r"[ ]*\[", "", value)
    if re.match(r"\][ ]*", value):
        value = re.sub(r"\][ ]*", "", value)
    return value


def remove_temp_and_output(value):
    if re.match(r"[ ]*temp\(", value):
        value = re.sub(r"[ ]*temp\(", "", value).rstrip(")")
    if re.match(r"[ ]*directory\(", value):
        value = re.sub(r"[ ]*directory\(", "", value).rstrip(")")
    return value


def markdown_table(rule_source, rule_schema):
    sections = ["input", "output", "params", "log", "benchmark",
                "resources", "threads", "conda", "container",
                "message", "shell", "wrapper", "script", "R"]

    def get_input_variabels(rule_source, include_section_regex, exclude_section_regex):
        rows = iter(rule_source.split("\n"))
        section_dict = {}
        for line in rows:
            if re.match(include_section_regex, line):
                line = next(rows)
                key = "missing"
                value = line
                if "=" in line:
                    key, value = re.split("[ ]*=[ ]*", line)
                    key = key.lstrip()
                for line in rows:
                    if re.match(exclude_section_regex, line):
                        section_dict[key] = replace_newline(remove_temp_and_output(remove_indent(value).rstrip(",")))
                        return section_dict
                    if "=" in line:
                        section_dict[key] = replace_newline(remove_temp_and_output(remove_indent(value).rstrip(",")))
                        key, value = re.split("[ ]*=[ ]*", line)
                        key = key.lstrip()
                    else:
                        value += line
                section_dict[key] = replace_newline(remove_temp_and_output(remove_indent(value).rstrip(",")))

        return section_dict

    header = "| Rule parameters | Key | Value | Description |" + \
             "\n| --- | --- | --- | --- |"
    table = ""

    def markdown_list(variables, variable_key, rule_schema, section_key, merge_section_sign):
        description = replace_newline(rule_schema['properties'][section_key]['properties'][variable_key]['description'])
        value = replace_newline(rule_schema['properties'][section_key]['properties'][variable_key].get('value', ""))
        # Variable is a list

        if value:
            variable_part_list = [value]
        else:
            if re.match(r"^[ ]\[", variables[variable_key]):
                variable_part_list = re.split(",[ ]*",
                                              re.sub(r"[, ]*\][ ]*$", "",
                                                     re.sub(r"^[ ]*\[", "",
                                                            variables[variable_key])))
            else:
                variable_part_list = [variables[variable_key]]
            variable_part_list = list(map(lambda v: f"`{v}`", variable_part_list))

        if merge_section_sign:
            section_key = merge_section_sign

        part_counter = len(variable_part_list)
        part_it = iter(variable_part_list)

        part = next(part_it)
        table = f"\n| {section_key} | {variable_key} | {part} | {description} |"

        if merge_section_sign is None and len(variables) > 1:
            section_key = "   "

        part_counter -= 1
        part_merge_sign = "  "
        if part_counter > 0:
            for part in part_it:
                part_counter -= 1
                if part_counter == 0:
                    part_merge_sign = "_ _"
                table += f"\n| {section_key} | {part_merge_sign} | {part} | {part_merge_sign} |"

        return table

    for section_key in rule_schema['properties']:
        section_copy = sections.copy()
        section_copy.remove(section_key)
        variables = get_input_variabels(rule_source,
                                        r"^[ \t]+" + section_key + ":",
                                        r"^[ \t]+(" + "|".join(section_copy) + r"):")
        count_variables = len(rule_schema['properties'][section_key]['properties'])
        variable_it = iter(rule_schema['properties'][section_key]['properties'])
        variable_key = next(variable_it)

        table += markdown_list(variables, variable_key, rule_schema, section_key, None)

        merge_section_sign = "   "
        if count_variables > 1:
            count_variables -= 1

            for variable_key in variable_it:
                count_variables -= 1

                if count_variables == 0:
                    merge_section_sign = "_  _"

                table += markdown_list(variables, variable_key, rule_schema, section_key, merge_section_sign)

    return header + table


def remove_trailing_empty_lines(rule):
    if rule.endswith("\n"):
        return remove_trailing_empty_lines(rule.rstrip())
    else:
        return rule


class markdown_gen:
    schema_file = ''
    tag = ''
    indent_val = "    "

    def safe_get_value(self, data, key):
        if data is None:
            return None, False
        try:
            output = data[key]
            return output, True
        except KeyError:
            return None, False

    def find_file(self, file_name):
        file_paths = list()

        if file_name is None:
            for folder_path in self.config_rule_folders:
                for file in glob.glob(f"{folder_path}/*.smk"):
                    file_paths.append(file)

            if file_paths is None:
                raise RuleFileMissingError(f"Could not locate any rule files at, {self.config_rule_folders}")
        else:
            for folder_path in self.config_rule_folders:
                file_path = os.path.join(folder_path, file_name)
                if os.path.exists(file_path):
                    file_paths.append(file_path)
                elif os.path.exists(f"{file_path}.smk"):
                    file_paths.append(f"{file_path}.smk")

            if file_paths is None:
                raise RuleFileMissingError(f"Could not locate rule file: {file_name}")

        return file_paths

    def extract_rule_source(self, file_name, rule_name):
        file_paths = self.find_file(file_name)

        for file_path in file_paths:
            if rule_name not in self.config_extracted_rules:
                rule_source = remove_trailing_empty_lines(extract_snakemake_rule(file_path, rule_name))
                if rule_source:
                    self.config_extracted_rules[rule_name] = {'source': rule_source, 'file_path': file_path}
                    return self.config_extracted_rules[rule_name]
            else:
                return self.config_extracted_rules[rule_name]
        if not self.config_extracted_rules[rule_name]:
            raise RuleMissingInSmkFileError("Not able to located {rule_name}")
        return None

    def add_schema_information(self, rule_name):
        if rule_name in self.config_extracted_rules and 'schema' in self.config_extracted_rules[rule_name]:
            return
        if rule_name not in self.config_extracted_rules:
            rule_source = self.extract_rule_source(None, rule_name)
            if not rule_source:
                raise RuleFileMissingError("table creation, unable to extrac source for rule {rule_name}")
            self.add_schema_information(rule_name)
        for schema in self.config_schemas:
            if rule_name in schema['properties']:
                self.config_extracted_rules[rule_name]['schema'] = schema['properties'][rule_name]
        if not self.config_extracted_rules[rule_name]['schema']:
            raise RuleMissingSchemaFileError(f"Rule name not found in provided schemas {self.config_schemas}")

    def get_markdown(self, markdown, **kwargs):
        for g in re.finditer(TAG_SNAKEMAKE_SOURCE_RE, markdown):
            parts = g.group()[1:-1].split("__")
            if len(parts) != 3:
                raise IncorrectTagError(f"Incorrect tag should be SNAKEMAKE_RULE_SOURCE__rulefilename__rulenamm was {g.group()}")
            _, file_name, rule_name = parts

            rule_source = self.extract_rule_source(file_name, rule_name)

            new_markdown = markdown_source(rule_source['source'])
            markdown = markdown.replace(g.group(), new_markdown)

        for g in re.finditer(TAG_SNAKEMAKE_TABLE_RE, markdown):
            parts = g.group()[1:-1].split("__")
            if len(parts) != 3:
                raise IncorrectTagError(f"Incorrect tag should be SNAKEMAKE_RULE_SOURCE__rulefilename__rulenamm was {g.group()}")
            _, file_name, rule_name = parts

            self.add_schema_information(rule_name)

            new_markdown = markdown_table(self.config_extracted_rules[rule_name]['source'],
                                          self.config_extracted_rules[rule_name]['schema'])

            markdown = markdown.replace(g.group(), new_markdown)

        return markdown

    def set_config(self, config):
        self.config_extracted_rules = dict()
        self.config_rule_folders = set()
        self.config_schemas = list()

        folders_config = self.safe_get_value(config, "rule_folders")
        if folders_config[0] and folders_config[1]:
            for folder in folders_config[0]:
                self.config_rule_folders.add(folder)

        schemas_config = self.safe_get_value(config, "schemas")
        if schemas_config[0] and schemas_config[1]:
            for schema in schemas_config[0]:
                with open(schema, 'r') as yaml_file:
                    self.config_schemas.append(yaml.safe_load(yaml_file))

        if not self.config_rule_folders:
            warnings.warn("No rule folders configured")
        if not self.config_schemas:
            warnings.warn("No rule schemas configured")
