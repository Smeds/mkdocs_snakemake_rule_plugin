import os
import re

TAG_SNAKEMAKE_SOURCE_RE = r"#SNAKEMAKE_RULE_SOURCE[A-Za-z0-9_-]+#"
TAG_SNAKEMAKE_TABLE_RE = r"#SNAKEMAKE_RULE_TABLE[A-Za-z0-9_-]+#"

class MissingConfigSettingError(ValueError):
    'Missing config setting'
    pass


class IncorrectTagError(ValueError):
    'Incorrect tag'
    pass


def extract_snakemake_rule_section(parts, data):
    if len(parts) > 0:
        return extract_snakemake_rule_section(parts[1:], data)
    else:
        return data

def remove_trailing_empty_lines(rule):
    if rule.endswith("\n"):
        return remove_trailing_empty_lines(rule.rstrip())
    else:
        return rule

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

    def get_markdown(self, markdown, **kwargs):
        for g in re.finditer(TAG_SNAKEMAKE_SOURCE_RE, markdown):
            parts = g.group()[1:-1].split("__")
            if len(parts) != 3:
                raise IncorrectTagError(f"Incorrect tag should be SNAKEMAKE_RULE_SOURCE__rulefilename__rulenamm was {g.group()}")            
            file_path = None
            for folder_path in self.config_rule_folders:
                print(folder_path)
                file_path = os.path.join(folder_path, parts[1])
                print(file_path)
                if os.path.exists(file_path):
                    break
                file_path += ".smk" 
                if os.path.exists(file_path):
                    break
                else:
                    file_path = None

            rule_source = remove_trailing_empty_lines(extract_snakemake_rule(file_path, parts[2]))
            new_markdown = self.markdown_source(rule_source)
            markdown = markdown.replace(g.group(), new_markdown)

        return markdown

    def markdown_source(self, code_item):
        return  "```\n" + code_item + "\n```" 
        
    def set_config(self, config):
        self.config_rule_folders = set()
        folders_config = self.safe_get_value(config, "rule_folders")
        if folders_config[0] and folders_config[1]:
            for folder in folders_config[0]:
                self.config_rule_folders.add(folder)
        if not self.config_rule_folders:
            raise MissingConfigSettingError("No rule folders configured")
