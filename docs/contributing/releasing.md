# Releasing

## How a release happens

1. Bump the version in `pkg/python/pyproject.toml` and `pkg/python/src/pymzlib/__init__.py`.
2. Update the changelog.
3. Tag and push:
   ```bash
   git tag v0.1.0 && git push origin v0.1.0
   ```
4. CI builds all four platform wheels, tests each on its real operating system, and — **if the
   `PYPI_PUBLISH` repository variable is `true`** — publishes to PyPI via
   [Trusted Publishing](https://docs.pypi.org/trusted-publishers/) — OIDC, so there are no API
   tokens stored anywhere.

!!! warning "PyPI publishing is off until you switch it on"
    The `publish to PyPI` job is gated on a repository variable `PYPI_PUBLISH=true` (Settings →
    Secrets and variables → Actions → Variables). It ships **disabled** on purpose, because PyPI
    cannot accept an upload until the two one-time steps below are done — a Trusted Publisher is
    registered and the file-size limit is granted. Until then, tagging safely builds and tests
    wheels on every platform without publishing, and a broken first version can't be shipped and
    locked. Turn it on only once both prerequisites are in place.

    **One-time PyPI setup, in order:**

    1. Register the Trusted Publisher: on PyPI, create/claim the `pymzlib` project, then add a
       GitHub Actions trusted publisher — owner `smith-chem-wisc`, repo `pyMzLib`, workflow
       `wheels.yml`, environment `pypi`. (See the file-size note below — do this together.)
    2. Request the file-size-limit increase (below); wait for it to be granted.
    3. Set the repository variable `PYPI_PUBLISH=true`.
    4. Bump to a real version and tag.

## Giving someone a build before the first PyPI release

You don't need PyPI to hand a colleague a working install. Every push to `main` builds all four
platform wheels; attach them to a GitHub Release and they install with one command:

```bash
# grab the wheels CI already built for a commit, then publish a release with them attached
gh run download <run-id> --pattern 'wheel-*' --dir wheels
gh release create v0.1.0.dev0 --prerelease --title "pyMzLib 0.1.0.dev0 (preview)" wheels/**/*.whl
```

The recipient installs the wheel for their platform straight from the release — no PyPI, no .NET:

```bash
pip install https://github.com/smith-chem-wisc/pyMzLib/releases/download/v0.1.0.dev0/<wheel-for-their-os>
```

Because `PYPI_PUBLISH` is off, this `v*` tag builds and tests wheels without attempting to publish.

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
