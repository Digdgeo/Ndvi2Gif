# Complete Jupyter Book Build Guide

This guide walks through building, testing, and deploying the Ndvi2Gif Jupyter Book tutorial.

## Prerequisites

### System Requirements

- Python 3.8+ (3.11 recommended)
- Git
- 2GB free disk space
- Internet connection (for downloading dependencies)

### Check Your Environment

```bash
python --version  # Should be 3.8+
git --version
```

## Step-by-Step Setup

### 1. Clone Repository (if not already done)

```bash
git clone https://github.com/Digdgeo/Ndvi2Gif.git
cd Ndvi2Gif
```

### 2. Create Virtual Environment (Recommended)

```bash
# Create environment
python -m venv book_env

# Activate
source book_env/bin/activate  # Linux/Mac
# OR
book_env\Scripts\activate  # Windows

# Verify
which python  # Should point to book_env
```

### 3. Install Dependencies

```bash
# Install Jupyter Book and dependencies
pip install -r book/requirements.txt

# Verify installation
jupyter-book --version
```

Expected output:
```
Jupyter Book      : 0.15.1
MyST-NB           : 0.17.2
Sphinx            : 5.3.0
```

## Building the Book

### First Build

```bash
# Clean any old builds
jupyter-book clean book/

# Build the book
jupyter-book build book/
```

You should see:
```
Running Jupyter-Book v0.15.1
Source Folder: book/
Config: book/_config.yml
Output Path: book/_build/html
...
===============================================================================

Finished generating HTML for book.
Your book's HTML pages are here:
    book/_build/html/
...
```

### Open in Browser

```bash
# Linux
xdg-open book/_build/html/index.html

# Mac
open book/_build/html/index.html

# Windows
start book/_build/html/index.html

# Or use Python HTTP server
cd book/_build/html
python -m http.server 8000
# Visit http://localhost:8000
```

## Development Workflow

### Making Changes

1. **Edit Content**
   ```bash
   # Edit a file
   nano book/tutorials/basic_ndvi.md
   # or use your favorite editor
   ```

2. **Rebuild**
   ```bash
   jupyter-book build book/
   ```

3. **Refresh Browser**
   - Reload the page to see changes

### Fast Iteration

For quick testing:

```bash
# Watch for changes (requires watchdog)
pip install watchdog
watchmedo shell-command \
    --patterns="*.md;*.ipynb" \
    --recursive \
    --command='jupyter-book build book/' \
    book/
```

### Check for Errors

Look for warnings/errors during build:

```bash
jupyter-book build book/ 2>&1 | grep -i "error\|warning"
```

Common issues:
- Missing files referenced in `_toc.yml`
- Broken internal links
- Invalid MyST syntax

## Testing Before Publishing

### 1. Check All Pages Load

```bash
# Start local server
cd book/_build/html
python -m http.server 8000
```

Visit http://localhost:8000 and click through all sections.

### 2. Validate Links

```bash
# Install linkchecker
pip install linkchecker

# Check links (this may take a while)
linkchecker --check-extern book/_build/html/index.html
```

### 3. Test Code Examples

Run code examples manually to verify they work:

```python
# From intro.md example
import ee
from ndvi2gif import NdviSeasonality

ee.Initialize()

roi = ee.Geometry.Point([-3.7, 40.4]).buffer(5000)
ndvi = NdviSeasonality(
    roi=roi,
    sat='S2',
    periods=12,
    start_year=2023,
    end_year=2024,
    index='ndvi'
)
print("✓ Example works!")
```

### 4. Check Build Size

```bash
du -sh book/_build/html
```

Should be < 100MB for efficient hosting.

## Publishing to GitHub Pages

### Option 1: Automatic Deployment (Recommended)

The `.github/workflows/deploy-book.yml` file automatically builds and deploys on push:

```bash
# Make changes
git add book/
git commit -m "Update Jupyter Book tutorials"
git push origin master

# GitHub Actions will build and deploy automatically
# Check progress at: https://github.com/Digdgeo/Ndvi2Gif/actions
```

### Option 2: Manual Deployment

```bash
# Build the book
jupyter-book clean book/
jupyter-book build book/

# Deploy to gh-pages branch
ghp-import -n -p -f book/_build/html

# Your book will be at:
# https://digdgeo.github.io/Ndvi2Gif/
```

### Setting Up GitHub Pages (First Time)

1. Go to repository settings
2. Navigate to "Pages" section
3. Source: Deploy from branch
4. Branch: `gh-pages`
5. Folder: `/ (root)`
6. Save

## Advanced Build Options

### Build as PDF

```bash
# Requires LaTeX
sudo apt-get install texlive-xetex texlive-fonts-recommended texlive-plain-generic  # Ubuntu/Debian

# Build PDF
jupyter-book build book/ --builder pdflatex

# Output: book/_build/latex/book.pdf
```

### Build Single Page

```bash
# Build only one page (faster for testing)
jupyter-book build book/ --path-output book/_build/singlepage tutorials/basic_ndvi.md
```

### Custom Builder Options

```bash
# Verbose output
jupyter-book build book/ -v

# Keep going despite errors
jupyter-book build book/ --keep-going

# Show all warnings
jupyter-book build book/ -W
```

## Configuration Tips

### Update Book Metadata

Edit `book/_config.yml`:

```yaml
title: Your Book Title
author: Your Name
logo: _static/logo.png  # Add your logo
```

### Add Custom CSS

```bash
# Create custom CSS
mkdir -p book/_static
cat > book/_static/custom.css << 'EOF'
/* Custom styles */
.admonition {
    border-radius: 5px;
}
EOF

# Reference in _config.yml
html:
  extra_css:
    - _static/custom.css
```

### Enable Features

In `_config.yml`:

```yaml
sphinx:
  config:
    # Add line numbers to code blocks
    highlight_options:
      linenostart: 1

    # Enable copy button
    copybutton_prompt_text: "$"
```

## Troubleshooting

### Issue: "Module not found" during build

```bash
pip install -r book/requirements.txt --upgrade
```

### Issue: Build hangs or is very slow

```bash
# Clean cache
jupyter-book clean book/ --all

# Disable notebook execution
# In _config.yml, ensure:
execute:
  execute_notebooks: off
```

### Issue: Pages not updating

```bash
# Force rebuild
jupyter-book clean book/ --all
jupyter-book build book/

# Clear browser cache (Ctrl+F5 or Cmd+Shift+R)
```

### Issue: Math equations not rendering

Check `_config.yml`:

```yaml
parse:
  myst_enable_extensions:
    - amsmath
    - dollarmath
```

### Issue: Images not showing

Use relative paths:

```markdown
# From tutorials/basic_ndvi.md to image in _static/
![Description](../_static/images/example.png)
```

### Issue: GitHub Actions deployment fails

Check:
1. GitHub Pages is enabled in repository settings
2. Workflow has write permissions
3. `book/requirements.txt` is complete
4. No broken links or missing files

## Best Practices

### 1. Test Before Committing

```bash
# Always test locally first
jupyter-book clean book/
jupyter-book build book/
# Check output for errors
```

### 2. Version Control

```bash
# Don't commit build artifacts
echo "_build/" >> book/.gitignore
echo ".jupyter_cache/" >> book/.gitignore
```

### 3. Keep It Fast

- Optimize images (< 500KB each)
- Set `execute_notebooks: off` for non-executable notebooks
- Avoid very large files

### 4. Maintain Structure

- Keep related content together
- Use meaningful file names
- Update `_toc.yml` when adding pages
- Follow naming conventions

## Maintenance

### Regular Updates

```bash
# Update dependencies quarterly
pip install --upgrade -r book/requirements.txt

# Test after updates
jupyter-book build book/
```

### Content Review

Periodically review:
- Broken links
- Outdated examples
- Code compatibility
- User feedback

## Getting Help

- **Jupyter Book Docs**: https://jupyterbook.org/
- **MyST Syntax**: https://myst-parser.readthedocs.io/
- **GitHub Issues**: https://github.com/executablebooks/jupyter-book/issues
- **Discussions**: https://github.com/executablebooks/meta/discussions

## Quick Reference

```bash
# Clean build
jupyter-book clean book/

# Build HTML
jupyter-book build book/

# Build and show warnings
jupyter-book build book/ -W

# Open in browser
xdg-open book/_build/html/index.html

# Deploy to GitHub Pages
ghp-import -n -p -f book/_build/html

# Check version
jupyter-book --version
```

---

Happy building! 🚀
