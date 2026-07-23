# Releasing

## How a release happens

1. Bump the version in `pkg/python/pyproject.toml` and `pkg/python/src/pymzlib/__init__.py`.
2. Update the changelog.
3. Tag and push:
   ```bash
   git tag v0.1.0 && git push origin v0.1.0
   ```
4. CI builds all four platform wheels, tests each on its real operating system, and publishes to
   PyPI via [Trusted Publishing](https://docs.pypi.org/trusted-publishers/) — OIDC, so there are
   no API tokens stored anywhere.

## Version numbers

`pymzlib` versions independently of mzLib. Tying them together would mean publishing a release
every time mzLib bumps a patch, which for a package on its 583rd version means a great deal of
noise for users whose API didn't change.

The mzLib commit a release was built from is recorded in `code/PINNED.md` and reported at runtime
by `pymzlib.bridge_version()`, so the provenance is always recoverable.

Semantic versioning, judged on the **Python** API: a change to the JSON envelope is internal
unless it changes what Python callers see.

## Before tagging

- [ ] `pytest` fully green, offline and network
- [ ] Wheel installs into a clean environment and runs (not the editable source tree)
- [ ] The `no-dotnet` CI job passes — the claim is untested without it
- [ ] `code/PINNED.md` matches the mzLib commit CI actually builds against
- [ ] Docs updated: coverage table, guide pages, changelog
- [ ] Bridge `protocol` bumped if the envelope changed incompatibly

## Channels

### PyPI — primary, automatic

Trusted Publishing on a tag. Nothing to do after setup.

!!! warning "The 100 MB file limit"
    PyPI rejects files over 100 MB by default, and pyMzLib's wheels are ~115 MB. Request an
    increase at [pypi.org/help](https://pypi.org/help/#file-size-limit) **before** the first
    release — it's routine and routinely granted (torch and friends all have one), but it isn't
    instant, and discovering it during a release is avoidable.

### bioconda — secondary, automatic after one-time setup

Bioconda is where proteomics tooling actually lives; pyOpenMS is there and not on conda-forge.

Submit a recipe once to [`bioconda-recipes`](https://github.com/bioconda/bioconda-recipes). After
it merges, bioconda's autobump bot watches PyPI and opens a pull request whenever a new version
appears — merge it and the package builds and publishes. As a bonus, every bioconda recipe
automatically gets a [BioContainer](https://biocontainers.pro/), so Docker distribution comes free.

!!! tip "Conda can do something pip cannot"
    conda-forge ships `dotnet-runtime` for linux-64, linux-aarch64, osx-64, osx-arm64 and win-64.
    A conda build could therefore declare `dotnet-runtime` as a dependency and ship a
    *framework-dependent* bridge — a package of a few megabytes instead of 115, with conda
    installing the runtime. The user experience is identical: one command, nothing to think about.

    That's a second build configuration and a second thing that can break, so it's worth doing
    only when bioconda is actually on the table — but it's the one packaging system where
    depending on .NET is clean rather than a burden pushed onto the user.

### Zenodo — for citation

Enable the [GitHub–Zenodo integration](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content)
before the first release. Every tagged release is then archived automatically and gets a DOI, and
the concept DOI cites "pyMzLib, any version" — which is what a methods section wants.

### conda-forge — probably not

Overlaps bioconda for this audience, and pyOpenMS's absence there is a hint. Each channel is a
recipe someone must maintain, and a stale recipe is somebody else's broken install.

## After a release

- Check the PyPI page renders (the README is the project description).
- Install from PyPI in a clean environment on at least one machine that isn't a CI runner.
- If the bioconda bot opened a PR, merge it.
