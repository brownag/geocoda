# Citation Guidelines for geocoda Documentation

## When to Add a Reference

Add a bibliography entry when:
1. A method/algorithm from a paper is implemented in geocoda code
2. Theoretical foundations require explanation beyond inline comments
3. Users would benefit from reading the original source for background
4. An R package is a core dependency (imports) or important optional tool (suggests)
5. External data sources (SSURGO, USDA standards) are used

Do NOT add references for:
- General R programming concepts
- Standard statistical methods without special compositional considerations
- Tangential topics not directly used by geocoda

## Bibliography Format (inst/references.bib)

All references must be in BibTeX format (.bib file). Use these templates based on reference type:

### Journal Article

```bibtex
@article{AuthorYear,
  author = {Last, F. M.},
  title = {Full Title in Title Case},
  journal = {Full Journal Name},
  volume = {X},
  number = {Y},
  pages = {start--end},
  year = {YYYY},
  doi = {10.xxxx/xxxxx},
  note = {[DETAILED ANNOTATION - see template below]}
}
```

### Book

```bibtex
@book{AuthorYear,
  author = {Last, F. M.},
  title = {Full Title in Title Case},
  publisher = {Publisher Name},
  year = {YYYY},
  address = {City, State/Country},
  edition = {Xth},  % if not 1st edition
  isbn = {XXX-X-XXX-XXXXX-X},
  doi = {10.xxxx/xxxxx},  % if available
  note = {[DETAILED ANNOTATION - see template below]}
}
```

### R Package

```bibtex
@manual{AuthorYear,
  title = {packagename: Package Title},
  author = {Last, F. M. and Last2, F. M.},
  year = {YYYY},
  note = {R package version X.Y.Z},
  url = {https://CRAN.R-project.org/package=packagename},
  note = {[DETAILED ANNOTATION - see template below]}
}
```

### Government/Institutional Resources

```bibtex
@misc{AuthorYear,
  author = {{Organization Name}},
  title = {Resource Title},
  year = {YYYY},
  publisher = {Publisher},
  url = {https://official.url/path}
}
```

## Annotation Requirements

All `note` fields MUST follow this structure:

```
note = {[WHAT IT PROVIDES] Theoretical/methodological contribution from this reference.
        [HOW GEOCODA USES IT] Specific geocoda functions, algorithms, or design decisions
        that rely on this work. Name functions explicitly (e.g., gc_ilr_params(),
        gc_sim_composition()).
        [USER BENEFIT] Why geocoda users should read this reference - what concepts will
        they understand better? What methodological choices will make sense?
        [ACCESS] Open access status (DOI link if OA), ISBN for books, or CRAN URL for
        R packages.}
```

**Guidelines:**
- **Length**: 3–6 sentences, 60–150 words per annotation
- **Specificity**: Include geocoda function names, chapter numbers, algorithm names
- **Clarity**: Explain WHY the reference matters to geocoda, not just WHAT it contains
- **Formatting**: Use `\%` for percent signs, `--` for dashes, proper capitalization

### Annotation Example

**Good annotation** (specific, actionable):
```bibtex
note = {Defines the Isometric Log-Ratio (ILR) transformation used throughout geocoda to
        convert constrained compositions (sand+silt+clay=100\%) to unconstrained coordinates
        suitable for geostatistics. Implemented in gc_ilr_params() for forward transformation
        and gc_sim_composition() for back-transformation via compositions::ilr(). Ensures
        sum constraints are automatically preserved after kriging/simulation. Essential for
        understanding why geocoda can guarantee valid soil textures. Open access via DOI.}
```

**Bad annotation** (too vague, non-specific):
```bibtex
note = {Important work on compositional data transformations}
```

## README and Vignette Citation Format

In Markdown (.Rmd) and plain Markdown (.md) files, use this format for references:

```markdown
- AuthorLastName, F. M. (Year). Title of Work. *Journal/Publisher*, Volume(Issue), pages.
  [https://doi.org/XX.XXXX/XXXXX](https://doi.org/XX.XXXX/XXXXX)
```

### Specific Rules

1. **Exact match**: Author names, year, title must match inst/references.bib entry exactly
2. **Include DOI**: Always add DOI hyperlink if available in bibliography (preferred over Journal URLs)
3. **Formatting**: Use `*Journal Name*` markdown for journal titles, `*Book Title*` for book titles
4. **Page ranges**: Use en-dash (`--`), not hyphen (`-`)
5. **Consistency**: All citations should use identical formatting throughout

### Examples

**Correct:**
```markdown
- Egozcue, J. J., et al. (2003). Isometric Log-Ratio Transformations for Compositional
  Data Analysis. *Mathematical Geology*, 35(3), 279–300.
  https://doi.org/10.1023/A:1023818214614
```

**Incorrect:**
```markdown
- Egozcue (2003) on ILR transformations
- Egozcue, J. J., Pawlowsky-Glahn, V., et al. Mathematical Geology 35: 279-300
```

## Vignette References Sections

Every vignette SHOULD include a References section with three subsections:

```markdown
## References

### Scientific Literature

[Papers and books providing theoretical foundations]
- Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
  Chapman and Hall, London. ISBN: 978-0-412-28060-3.

### R Packages

[Core packages used in examples or by geocoda]
- Pebesma, E. J. (2004). Multivariable geostatistics in S: the gstat package.
  *Computers & Geosciences*, 30, 683–691. https://doi.org/10.1016/j.cageo.2004.03.012

### Function Documentation

[For convenience, list key geocoda functions documented in vignette]
- `?gc_ilr_params` - Estimate ILR transformation parameters
- `?gc_sim_composition` - Generate constrained spatial realizations
```

## URL Guidelines

### Best Practices

- **DOIs** (preferred): Always use DOI links (https://doi.org/XXXX/XXXXX) - they are permanent and always resolve
- **R packages**: Use standard format `https://CRAN.R-project.org/package=name` for CRAN packages
- **Government**: Use official .gov URLs for USDA/NRCS resources (not mirrors or caches)
- **GitHub**: Use `https://github.com/user/repo` for development packages
- **Journal URLs**: Provide as fallback if DOI unavailable
- **Author preprints**: Do NOT include ResearchGate, arXiv, or personal preprint URLs

### Specific Resources

**CRAN packages:**
```
https://CRAN.R-project.org/package=packagename
```

**USDA Soil Survey Manual:**
```
https://www.nrcs.usda.gov/resources/guides-and-instructions/soil-survey-manual
```

**SSURGO Database:**
```
https://www.nrcs.usda.gov/resources/data-and-reports/soil-survey-geographic-ssurgo-database
```

**Soil Data Access API:**
```
https://sdmdataaccess.nrcs.usda.gov/
```

## Verification Checklist

Before committing reference changes, ensure:

- [ ] All entry types (article, book, manual, etc.) are correct
- [ ] All required fields present (author, title, year, journal/publisher)
- [ ] All DOIs resolve correctly (test with `curl -I https://doi.org/[DOI]`)
- [ ] All URLs return HTTP 200 (test with `curl -I [URL]`)
- [ ] Annotations are 3–6 sentences with specific geocoda function names
- [ ] README.Rmd citations match inst/references.bib entries exactly
- [ ] All BibTeX syntax is valid (check with bibtex::read.bib())
- [ ] No spelling errors in titles or author names
- [ ] ISBN numbers are correctly formatted

## Annual Maintenance

At least once per year (ideally before each major release):

1. **Test all DOIs**: Verify DOI links still resolve correctly
2. **Check R package versions**: Update year and version numbers for package citations
3. **Verify government URLs**: Confirm USDA/NRCS links still point to correct resources
4. **Check for broken links**: Run automated link checker on README and vignettes
5. **Update deprecations**: Remove references to deprecated packages or methods

## For Maintainers

### Adding a New Reference

1. Check if reference already exists in inst/references.bib
2. Determine reference type (article, book, manual, misc)
3. Gather all bibliographic information (DOI preferred over URL)
4. Write detailed annotation following template above
5. Add entry to inst/references.bib in alphabetical order by citation key
6. Update README.Rmd or relevant vignette References section
7. Test with `bibtex::read.bib("inst/references.bib")`
8. Verify DOI and URLs resolve correctly
9. Regenerate README.md if README.Rmd was edited

### Example: Adding a New Article

```bibtex
@article{NewAuthor2024,
  author = {New, First M. and Author, Second A.},
  title = {Compositional Geostatistics in Practice},
  journal = {Journal of Applied Geostatistics},
  volume = {15},
  number = {2},
  pages = {123--145},
  year = {2024},
  doi = {10.1234/jag.2024.001},
  note = {Advanced practical guide for applying compositional kriging methods.
          Demonstrates methods used in gc_fit_vgm() for ILR variogram modeling
          and gc_sim_composition() for ensemble simulation. Includes case studies
          on soil texture mapping. Helpful for practitioners implementing geocoda
          workflows on large spatial datasets. Open access via DOI.}
}
```

## Citation Key Naming Convention

Use the format: `AuthorLastName` + `Year` (or abbreviated if multiple authors)

Examples:
- `Aitchison1986`
- `Egozcue2003`
- `PalareaaAlbaladejo2015`  (use first author surname)
- `USDA2017` (for government documents)
- `NRCSSSURGO` (for major databases/systems)

## Questions or Issues?

If you have questions about citations or references:
1. Check this document first
2. Look at existing entries in inst/references.bib for examples
3. Review recent vignette References sections for current standards
4. Open an issue on GitHub if guidance is needed
