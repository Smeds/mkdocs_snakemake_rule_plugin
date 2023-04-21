from mkdocs import plugins, config
from .markdown import markdown_gen


class SnakemakeRule(plugins.BasePlugin):

    config_scheme = (
        ("rule_folders", config.config_options.Type(list, default=[])),
        ("schemas", config.config_options.Type(list, default=[])),
    )

    generator = markdown_gen()

    def on_config(self, config):
        self.generator = markdown_gen()
        self.generator.set_config(self.config)

    def __init__(self):
        super().__init__()

    def on_page_markdown(self, markdown, **kwargs):
        return self.generator.get_markdown(markdown, **kwargs)
