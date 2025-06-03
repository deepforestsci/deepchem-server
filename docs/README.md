# DeepChem Server Documentation

This directory contains the complete documentation for the DeepChem Server project, built with [Sphinx](https://www.sphinx-doc.org/) and formatted in [reStructuredText (RST)](https://docutils.sourceforge.io/rst.html).

## Documentation Structure

```
docs/
├── Makefile                    # Build automation
├── README.md                   # This file
├── source/                     # Documentation source files
│   ├── conf.py                # Sphinx configuration
│   ├── index.rst              # Main documentation index
│   ├── get_started/           # Getting started guide
│   │   ├── index.rst
│   │   ├── installation.rst
│   │   ├── quick_start.rst
│   │   ├── configuration.rst
│   │   └── examples.rst
│   ├── api_reference/         # API documentation
│   │   ├── index.rst
│   │   ├── rest_api.rst
│   │   ├── core_modules.rst
│   │   ├── data_structures.rst
│   │   └── utils.rst
│   └── py_ds_library/         # Python client library docs
│       └── index.rst
└── build/                     # Generated documentation (ignored by git)
```

## Quick Start

### Prerequisites

1. **Python 3.8+** with pip
2. **Sphinx and dependencies**:
   ```bash
   pip install sphinx sphinx-rtd-theme
   ```

### Building Documentation

1. **Navigate to the docs directory**:
   ```bash
   cd docs
   ```

2. **Build HTML documentation**:
   ```bash
   make html
   ```

3. **View the documentation**:
   ```bash
   make serve
   ```
   Then open http://localhost:8000 in your browser.

## Available Make Targets

### Basic Targets
- `make html` - Build HTML documentation
- `make clean` - Remove all build artifacts
- `make help` - Show all available targets

### Development Targets
- `make serve` - Build and serve documentation locally on port 8000
- `make watch` - Auto-rebuild documentation when files change (requires `sphinx-autobuild`)
- `make fast` - Fast build (skips some expensive operations)
- `make strict` - Build with warnings treated as errors
- `make rebuild` - Clean build directory and rebuild everything

### Quality Assurance
- `make linkcheck` - Check for broken external links
- `make check-internal` - Check for broken internal links only
- `make coverage` - Generate documentation coverage report
- `make spelling` - Check spelling (requires `sphinxcontrib-spelling`)

### Other Formats
- `make latex` - Generate LaTeX output
- `make latexpdf` - Generate PDF via LaTeX
- `make epub` - Generate EPUB e-book
- `make man` - Generate manual pages
- `make text` - Generate plain text output

### Utility Targets
- `make stats` - Show documentation statistics
- `make apidoc` - Generate API documentation from source code
- `make dev-setup` - Show development setup instructions

## GitHub Workflows

The documentation includes automated GitHub workflows for testing:

### Documentation Build Testing (`docs.yml`)

Automatically runs on:
- Push to `main` or `develop` branches
- Pull requests to `main` or `develop` branches
- Manual trigger via GitHub Actions UI

**Features:**
- **Build Testing**: Tests documentation builds in both strict and permissive modes
- **Style Checking**: Lints RST files with doc8 for style consistency
- **Makefile Testing**: Validates all Makefile targets work correctly
- **Dependency Caching**: Speeds up builds by caching pip dependencies

**Jobs:**
1. `docs-build` - Main documentation build with strict error checking
2. `test-makefile` - Validates Makefile targets work correctly

### Read the Docs Deployment

Documentation is automatically deployed to Read the Docs when changes are pushed to the main branch. Read the Docs will:

- **Automatically build** documentation from the latest main branch
- **Handle hosting** with professional documentation hosting features
- **Provide search functionality** across all documentation
- **Support multiple versions** if configured
- **Offer PDF/ePub downloads** automatically

**Read the Docs Benefits:**
- Professional documentation hosting optimized for technical docs
- Automatic builds from your Git repository
- Built-in search functionality
- Version management (latest, stable, feature branches)
- Download options (PDF, ePub, HTML zip)
- Analytics and usage statistics
- Custom domain support

**Configuration:**
Read the Docs automatically detects Sphinx projects and uses the `docs/` directory structure. The `docs/requirements.txt` file ensures all dependencies are installed during the build process.

### Viewing Documentation

- **Live Documentation**: Available at your Read the Docs URL
- **Build Status**: Check the Actions tab in your GitHub repository for testing
- **Documentation Artifacts**: Download built HTML from workflow runs for preview
- **Build Logs**: Review detailed logs for any build failures

### Local Testing Before Push

Before pushing documentation changes, test locally to ensure they'll build successfully on Read the Docs:

```bash
# Test the same build process as CI and Read the Docs
cd docs
make clean
make html

# Test Makefile targets
make help
make stats

# Check for style issues (install doc8 first)
pip install doc8
doc8 source/ --max-line-length=120 --allow-long-titles
```

**Note:** The GitHub workflow tests ensure your documentation builds correctly before it reaches Read the Docs, preventing broken builds from being deployed.

## Development Workflow

### Setting Up Development Environment

1. **Install required packages**:
   ```bash
   pip install sphinx sphinx-rtd-theme sphinx-autobuild
   ```

2. **Optional packages for enhanced features**:
   ```bash
   pip install sphinxcontrib-spelling
   ```

3. **Check your setup**:
   ```bash
   make dev-setup
   ```

### Making Changes

1. **Start auto-rebuild server** (recommended for development):
   ```bash
   make watch
   ```
   This will automatically rebuild and refresh your browser when you make changes.

2. **Edit documentation files** in the `source/` directory using reStructuredText format.

3. **Preview changes** by visiting http://localhost:8000

### Adding New Documentation

1. **Create new `.rst` files** in the appropriate directory under `source/`
2. **Add references** to new files in the relevant `index.rst` or other files
3. **Update navigation** by modifying the `toctree` directives as needed

### Code Documentation

The documentation automatically includes API documentation generated from Python docstrings in the source code. To update:

1. **Add/update docstrings** in the source code using Google or NumPy style
2. **Regenerate API docs**:
   ```bash
   make apidoc
   ```
3. **Rebuild documentation**:
   ```bash
   make html
   ```

## Documentation Conventions

### File Naming
- Use lowercase with underscores: `quick_start.rst`
- Use descriptive names that reflect content
- Keep filenames concise but clear

### reStructuredText Style
- Use consistent heading hierarchy:
  ```rst
  Main Title
  ==========
  
  Section
  -------
  
  Subsection
  ^^^^^^^^^^
  
  Sub-subsection
  """"""""""""""
  ```

### Code Examples
- Include working code examples with proper syntax highlighting
- Test code examples to ensure they work
- Use realistic examples that users can run

### Cross-References
- Use Sphinx cross-references for internal links:
  ```rst
  :doc:`installation`
  :ref:`section-label`
  :class:`ClassName`
  :func:`function_name`
  ```

## Troubleshooting

### Common Issues

1. **Build fails with import errors**:
   - Ensure the DeepChem Server dependencies are installed
   - Check that the source code paths in `conf.py` are correct

2. **Links are broken**:
   - Run `make linkcheck` to identify broken links
   - Update or remove broken external links

3. **Documentation looks wrong**:
   - Clear the build cache: `make clean-cache`
   - Rebuild: `make rebuild`

4. **Auto-rebuild not working**:
   - Install sphinx-autobuild: `pip install sphinx-autobuild`
   - Check that you're in the `docs/` directory when running `make watch`

### Getting Help

- **Sphinx Documentation**: https://www.sphinx-doc.org/
- **reStructuredText Guide**: https://docutils.sourceforge.io/rst.html
- **Read the Docs Theme**: https://sphinx-rtd-theme.readthedocs.io/

## Contributing

When contributing to the documentation:

1. Follow the existing structure and conventions
2. Test your changes locally before submitting
3. Run quality checks: `make linkcheck` and `make strict`
4. Update this README if you add new features or change workflows

## Publishing

The documentation can be deployed to various platforms:

- **Read the Docs**: Connect your repository to automatically build and host
- **GitHub Pages**: Use the `make github` target (if configured)
- **Static hosting**: Upload the contents of `build/html/` to any web server

For automated deployment, consider setting up CI/CD pipelines that run:
```bash
cd docs
make clean
make strict  # Build with warnings as errors
make linkcheck  # Verify links
``` 