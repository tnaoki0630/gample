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

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', "sphinx.ext.viewcode", "myst_parser"]

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
    pattern = (r'-\s*\(NSMutableDictionary\*\)\s*(\w+)\s*\{.*?'
               r'return\s*\[\s*@\{\s*(.*?)\s*\}\s*mutableCopy\s*\];')
    blocks = re.findall(pattern, text, flags=re.S)

    value_pat = r'(?:@\[[^\]]*\]|@\([^)]*\)|@"[^"]*"|@[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?|\w+)'
    lines = ["defaults", "==============", ""]
    for name, body in blocks:
        items = re.findall(r'@"([^"]+)"\s*:\s*(' + value_pat + r')', body)
        lines += [name, "-"*len(name), "",
                  ".. list-table::", "   :header-rows: 1", "",
                  "   * - Key", "     - Default"]
        for k, v in items:
            lines += [f"   * - ``{k}``", f"     - ``{v.strip()}``"]
        lines += [""]

    (OUT / "defaults.rst").write_text("\n".join(lines), encoding="utf-8")
    print(f"blocks_found={len(blocks)}")

def setup(app):
    app.connect("builder-inited", build_defaults_page)

# -- Auto generation of bibliography --------
extensions += ["sphinxcontrib.bibtex"]
# 単一ファイルなら:
# bibtex_bibfiles = ["bibliography/refs.bib"]
# 複数/ワイルドカードで読み込みたい場合:
from glob import glob
bibtex_bibfiles = sorted(glob("bibliography/*.bib"))
bibtex_default_style = "unsrt"          # 例: 並びは出現順
bibtex_reference_style = "author_year"  # 例: (著者, 年)
print(f"bibtex: {bibtex_bibfiles}")

# -- mermaid w/ TeX setting ----------------
extensions += ['sphinxcontrib.mermaid']
# Mermaid 10.9+ を指定（数式サポートのため）
mermaid_version = "11.12.0"
# 必要なら KaTeX に強制（各環境で表示差を避ける）
mermaid_init_js = "mermaid.initialize({startOnLoad:true, forceLegacyMathML:true});"
# KaTeX の CSS を供給（forceLegacyMathML 使用時）
html_css_files = [
  "https://cdn.jsdelivr.net/npm/katex@0.16.10/dist/katex.min.css",
]