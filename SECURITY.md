# Security policy

## Reporting a vulnerability

Please report security issues **privately**, not as a public issue:

- [Open a private security advisory](https://github.com/smith-chem-wisc/pyMzLib/security/advisories/new)
  (preferred), or
- email **mm_support@chem.wisc.edu**

Please include what you found, how to reproduce it, and what an attacker could achieve. We'll
acknowledge within a week and keep you informed while we work on it. If you'd like credit in the
advisory, say so and tell us how you'd like to be named.

## What's in scope

pyMzLib runs a bundled executable and fetches data from public repositories, so the interesting
surface is:

- **The bundled .NET payload** — a vulnerability in the runtime or in mzLib's dependencies that
  reaches users through our wheel. We rebuild against a patched runtime and release.
- **Path handling on download.** Filenames come from PRIDE and are treated as untrusted; the
  bridge refuses any name that isn't a bare leaf name, so a crafted filename cannot write outside
  the destination directory. A way around that check is a vulnerability — please report it.
- **The bridge's argument handling**, if a way exists to make pyMzLib execute something other
  than its own bundled payload.
- **Release integrity** — anything that would let a third party publish a package as us.

## What's not in scope

- Vulnerabilities in mzLib itself. Please report those to
  [mzLib](https://github.com/smith-chem-wisc/mzLib/security) — though if you're unsure which side
  a problem sits on, report it here and we'll route it.
- Denial of service against EBI's PRIDE Archive, which isn't ours.
- The fact that `download()` writes files you asked it to write.

## Supported versions

pyMzLib is pre-1.0 and moves fast. Fixes go into the latest release; there are no long-term
support branches yet. Once 1.0 ships, this section will say something more useful.

## A note on the bundled runtime

Because each wheel contains a complete .NET runtime, a runtime security patch reaches you through
a pyMzLib release rather than through your system package manager. That's a real obligation on
us: we track .NET servicing releases and rebuild. If you notice we're behind on one, that's worth
an issue — a normal public one is fine.
