# Contributing to DiplotypeForge

DiplotypeForge is developed by the Society for Immune Polymorphism (SIP) using
the NMDP Bioinformatics **fork-and-pull-request** workflow. Contributions are
made from a personal fork and merged into the canonical SIP repository through
reviewed pull requests.

- Canonical upstream: <https://github.com/societyforimmunepolymorphism/DiplotypeForge>
- Reference workflow (NMDP `dash`):
  <https://github.com/nmdp-bioinformatics/dash/blob/master/README.git.fork.method.md>

> **Adaptation note:** the NMDP reference document uses `master` as the default
> branch. DiplotypeForge uses **`main`**; substitute `main` for `master`
> throughout.

## Fork method (adapted to `main`)

1. **Fork** `societyforimmunepolymorphism/DiplotypeForge` into your personal
   GitHub account using the GitHub *Fork* button (this records the parent/fork
   relationship).
2. **Clone your personal fork** as `origin`:
   ```bash
   git clone https://github.com/<your-username>/DiplotypeForge.git
   cd DiplotypeForge
   ```
3. **Configure SIP as `upstream`:**
   ```bash
   git remote add upstream https://github.com/societyforimmunepolymorphism/DiplotypeForge.git
   ```
4. **Synchronize local `main` with `upstream/main`:**
   ```bash
   git checkout main
   git pull upstream main
   ```
5. **Push the synchronized `main` to your fork:**
   ```bash
   git push origin main
   ```
6. **Create a feature branch:**
   ```bash
   git checkout -b my-feature
   ```
7. **Make and test changes** (see *Development* below).
8. **Commit and push the feature branch to your fork:**
   ```bash
   git pull upstream main          # pick up any upstream changes first
   git commit -m "describe the change"
   git push origin my-feature
   ```
9. **Open a pull request** against `societyforimmunepolymorphism/DiplotypeForge`
   `main` from your fork's feature branch.
10. **Address review feedback** by pushing additional commits to the same
    feature branch.
11. **After merge**, delete the now-obsolete remote and local feature branches
    and resynchronize `main`:
    ```bash
    git push origin --delete my-feature
    git checkout main && git pull upstream main && git push origin main
    git branch -d my-feature
    ```

Never force-push shared branches, and never rewrite published history.

## Development

Requires Python 3.12+ and [uv](https://docs.astral.sh/uv/).

```bash
uv venv --python 3.12
uv pip install -e ".[dev]"      # add ,viz for the matplotlib experiment plot
```

Before opening a pull request, ensure:

```bash
uv run ruff check src tests examples      # lint — must be clean
uv run pytest -q                          # tests — must pass
uv build                                  # sdist + wheel build must succeed
```

- Match the existing code style and conventions; keep changes scoped.
- Add tests for new behavior and bug fixes.
- Preserve reproducibility guarantees (seed hierarchy, deterministic
  serialization). The CI workflow runs ruff, pytest, a build, a clean-environment
  wheel install, and a representative CLI simulation on Python 3.12 and 3.13.
- Do not change the project license without explicit SIP approval.

## Scope

DiplotypeForge simulates known, phased diploid genotypes from population
haplotype-frequency distributions. Laboratory typing-observation simulation is
out of scope and belongs to the companion project **TissueTypeForge** (see
`docs/DIPLOTYPEFORGE_AND_TISSUETYPEFORGE.md`). The supported public API is
defined in `docs/API_STABILITY.md`.
