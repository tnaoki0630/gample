# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'gample'
language = "ja"
copyright = '2025, Naoki Tsunezawa'
author = 'Naoki Tsunezawa'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', "sphinx.ext.viewcode"]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Auto generation of solver defaults -------------------------------------------------
from pathlib import Path
import re

SRC = Path(__file__).parents[2] / "src" / "Init.mm"
OUT = Path(__file__).parent / "_autogen"
OUT.mkdir(exist_ok=True)

def build_defaults_page(app):
    text = SRC.read_text(encoding="utf-8")
    blocks = re.findall(
        r'-\s*\(NSMutableDictionary\*\)\s*(\w+)\s*\{.*?@(\{.*?\})\s*\}\s*mutableCopy\]',
        text, flags=re.S)
    lines = ["デフォルト値一覧", "==============", ""]
    for name, body in blocks:
        # キー:値 を素朴に抽出（@"K": @V / @"K": @"str" / @"K": @[...])
        items = re.findall(r'@\"([^\"]+)\"\s*:\s*([^,}]+)', body)
        lines += [f"{name}", "-"*len(name), ""]
        lines += [".. list-table::", "   :header-rows: 1", "", "   * - Key", "     - Default",]
        for k, v in items:
            lines += [f"   * - ``{k}``", f"     - ``{v.strip()}``"]
        lines += [""]
    (OUT / "defaults.rst").write_text("\n".join(lines), encoding="utf-8")

def setup(app):
    app.connect("builder-inited", build_defaults_page)

# -- Auto generation of bibliography --------
extensions += ["sphinxcontrib.bibtex", "myst_parser"]
# 単一ファイルなら:
# bibtex_bibfiles = ["bibliography/refs.bib"]
# 複数/ワイルドカードで読み込みたい場合:
from glob import glob
bibtex_bibfiles = sorted(glob("bibliography/*.bib"))
bibtex_default_style = "unsrt"          # 例: 並びは出現順
bibtex_reference_style = "author_year"  # 例: (著者, 年)