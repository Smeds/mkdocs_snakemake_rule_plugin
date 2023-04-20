
# -*- coding: utf-8 -*-
from setuptools import setup
import versioneer

with open("README.md") as f:
    readme = f.read()

setup(
    name="mkdocs-snakemake-rule-plugin",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="MkDocs Plugin used to extrace rule information, to generate code and tables.",
    long_description_content_type="text/markdown",
    long_description=readme,
    keywords=["mkdocs", "plugin", "yaml", "schema"],
    author="Patrik Smeds",
    author_email="patrik.smeds@gmail.com",
    url="https://github.com/smeds/mkdocs-snakemake-rule-plugin",
    license="MIT license",
    packages=["mkdocs_snakemake_rule_plugin"],
    install_requires=["mkdocs", "pyyaml"],
    python_requires=">=3.8, <4",
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    # This entry point is necessary for MkDocs to be able to use the plugin
    entry_points={
        'mkdocs.plugins': [
            'snakemake-rule = mkdocs_snakemake_rule_plugin.plugin:SnakemakeRule',
        ]
    },
)
