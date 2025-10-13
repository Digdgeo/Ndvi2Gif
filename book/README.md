# Ndvi2Gif Jupyter Book

This directory contains the source files for the Ndvi2Gif tutorial Jupyter Book.

## Building the Book Locally

### Prerequisites

```bash
# Install Jupyter Book and dependencies
pip install -r requirements.txt
```

### Build HTML

```bash
# Clean previous build
jupyter-book clean book/

# Build the book
jupyter-book build book/
```

The HTML output will be in `book/_build/html/`. Open `book/_build/html/index.html` in your browser.

### Build PDF (Optional)

```bash
jupyter-book build book/ --builder pdflatex
```

Requires LaTeX installation.

## Directory Structure

```
book/
├── _config.yml              # Jupyter Book configuration
├── _toc.yml                 # Table of contents
├── intro.md                 # Landing page
├── getting_started/         # Installation and setup
├── tutorials/               # Step-by-step tutorials
├── advanced/                # Advanced features
├── use_cases/               # Real-world examples
├── notebooks/               # Jupyter notebooks
├── reference/               # API docs and reference
└── _build/                  # Generated output (gitignored)
```

## Publishing to GitHub Pages

The book is automatically built and published via GitHub Actions on every push to `master`.

Manual publishing:

```bash
# Build the book
jupyter-book build book/

# Publish to gh-pages branch
ghp-import -n -p -f book/_build/html
```

## Contributing to Documentation

1. Edit markdown files in appropriate directories
2. For new pages, add to `_toc.yml`
3. Build locally to test: `jupyter-book build book/`
4. Submit a pull request

### Writing Style Guidelines

- Use clear, concise language
- Include code examples for all features
- Add tips, notes, and warnings using MyST syntax
- Test all code examples
- Include visualizations where helpful

### MyST Markdown Syntax

```markdown
# Admonitions
:::{note}
This is a note
:::

:::{tip}
This is a tip
:::

:::{warning}
This is a warning
:::

# Code blocks with syntax highlighting
\`\`\`python
import ndvi2gif
\`\`\`

# Math
$$
NDVI = \frac{NIR - Red}{NIR + Red}
$$

# Cross-references
[Link text](path/to/page.md)

# Images
\`\`\`{image} path/to/image.png
:alt: Alt text
:width: 600px
\`\`\`
```

## Updating Notebooks

Notebooks in `notebooks/` should:
1. Be fully executed with outputs
2. Include clear markdown explanations
3. Use reproducible examples
4. Follow PEP 8 style
5. Include error handling examples

## Troubleshooting Build Issues

### Issue: Missing dependencies

```bash
pip install -r requirements.txt --upgrade
```

### Issue: Kernel errors in notebooks

Set `execute_notebooks: off` in `_config.yml` to skip execution.

### Issue: Math rendering

Ensure `amsmath` extension is enabled in `_config.yml`.

## Resources

- [Jupyter Book Documentation](https://jupyterbook.org/)
- [MyST Markdown Syntax](https://myst-parser.readthedocs.io/)
- [Sphinx Documentation](https://www.sphinx-doc.org/)
