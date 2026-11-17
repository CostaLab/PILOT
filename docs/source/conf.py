# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'PILOT'
copyright = '2023, Costalab'
author = 'Mehdi Joodaki & Mina Shaigan'

release = '2.1.1'
version = '2.1.1'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'myst_parser',
    'nbsphinx',  # for jupyter notebooks
    'sphinx_design',  # cards and grids
]

myst_enable_extensions = ['colon_fence']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']
exclude_patterns = ['_build']

# -- Options for HTML output

html_theme = 'sphinx_book_theme'
html_title = 'PILOT'
html_logo = 'logo.png'
html_favicon = 'logo.png'
html_static_path = ['_static']
html_css_files = ['custom.css']

html_theme_options = {
    'repository_url': 'https://github.com/CostaLab/PILOT',
    'repository_branch': 'main',
    'path_to_docs': 'docs/source',
    'use_repository_button': True,
    'use_issues_button': True,
    'use_download_button': True,
    'home_page_in_toc': True,
    'show_toc_level': 2,
    'navigation_with_keys': False,
    'logo': {'text': f'PILOT {release}'},
    'icon_links': [
        {'name': 'PyPI', 'url': 'https://pypi.org/project/pilotpy/',
         'icon': 'fa-brands fa-python', 'type': 'fontawesome'},
    ],
}

# -- Options for EPUB output
epub_show_urls = 'footnote'
