# Jupyter Book Quick Start

Get your Jupyter Book tutorial up and running in minutes!

## Installation

```bash
cd /home/diego/git/Ndvi2Gif
pip install -r book/requirements.txt
```

## Build the Book

```bash
# Clean any previous builds
jupyter-book clean book/

# Build the book
jupyter-book build book/

# Open in browser
firefox book/_build/html/index.html
# or
google-chrome book/_build/html/index.html
```

## Development Workflow

### 1. Edit Content

Edit markdown files in `book/`:
- `intro.md` - Landing page
- `getting_started/*.md` - Installation and setup guides
- `tutorials/*.md` - Tutorial content
- `advanced/*.md` - Advanced features
- `reference/*.md` - API documentation

### 2. Add New Pages

To add a new page:

1. Create the markdown file:
   ```bash
   touch book/tutorials/my_new_tutorial.md
   ```

2. Add to table of contents in `book/_toc.yml`:
   ```yaml
   - caption: Core Tutorials
     chapters:
       - file: tutorials/my_new_tutorial
         title: "My New Tutorial"
   ```

3. Rebuild the book

### 3. Test Locally

Always test before pushing:

```bash
jupyter-book build book/
```

Check for errors in the output.

### 4. Preview in Browser

```bash
# After building
cd book/_build/html
python -m http.server 8000
# Visit http://localhost:8000
```

## Publishing

### Automatic (Recommended)

The book automatically builds and deploys to GitHub Pages when you push to the `master` branch:

```bash
git add book/
git commit -m "Update Jupyter Book content"
git push origin master
```

Visit: `https://digdgeo.github.io/Ndvi2Gif/`

### Manual Publishing

```bash
# Build
jupyter-book build book/

# Publish to gh-pages
ghp-import -n -p -f book/_build/html
```

## Common Tasks

### Add a New Tutorial

```bash
# Create file
cat > book/tutorials/my_tutorial.md << 'EOF'
# My Tutorial Title

Introduction...

## Section 1

Content...
EOF

# Add to _toc.yml
# Then rebuild
jupyter-book build book/
```

### Add Images

```bash
# Create images directory
mkdir -p book/_static/images

# Add image to markdown
# ![Alt text](../_static/images/my_image.png)
```

### Add Math

Use LaTeX syntax:

```markdown
Inline math: $E = mc^2$

Display math:
$$
NDVI = \frac{NIR - Red}{NIR + Red}
$$
```

### Add Code Blocks

```markdown
\`\`\`python
import ndvi2gif
print("Hello, World!")
\`\`\`
```

### Add Admonitions

```markdown
:::{note}
This is a note
:::

:::{tip}
This is a helpful tip
:::

:::{warning}
This is a warning
:::
```

## Directory Structure

```
book/
├── _config.yml              # Main configuration
├── _toc.yml                 # Table of contents
├── intro.md                 # Landing page
├── getting_started/         # Installation guides
│   ├── installation.md
│   ├── authentication.md
│   ├── quick_start.md
│   └── roi_options.md
├── tutorials/               # Core tutorials
│   ├── basic_ndvi.md
│   ├── multi_sensor.md
│   └── ...
├── advanced/                # Advanced features
├── use_cases/               # Real-world examples
├── notebooks/               # Jupyter notebooks
├── reference/               # API docs
│   ├── api.md
│   ├── indices.md
│   ├── faq.md
│   └── references.bib
├── _static/                 # Static assets (images, CSS)
└── _build/                  # Generated output (gitignored)
```

## Troubleshooting

### Build fails with "missing module"

```bash
pip install -r book/requirements.txt --upgrade
```

### Links are broken

Check that:
- File paths are correct in `_toc.yml`
- Cross-references use correct relative paths
- All referenced files exist

### Math not rendering

Ensure `amsmath` is enabled in `_config.yml`:

```yaml
parse:
  myst_enable_extensions:
    - amsmath
    - dollarmath
```

### Images not showing

Use relative paths from the markdown file:

```markdown
# From tutorials/basic_ndvi.md
![Image](../_static/images/example.png)
```

## Next Steps

1. **Customize**: Edit `_config.yml` to update title, author, URLs
2. **Write Content**: Start with high-priority tutorials
3. **Add Notebooks**: Convert existing notebooks in `examples_notebooks/`
4. **Review**: Test all links and code examples
5. **Deploy**: Push to GitHub and enable Pages

## Resources

- [Jupyter Book Docs](https://jupyterbook.org/)
- [MyST Markdown Guide](https://myst-parser.readthedocs.io/)
- [Example Books](https://executablebooks.org/en/latest/gallery.html)

---

Happy documenting! 📚
