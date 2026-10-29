---
name: land
description: >-
  Land explicitly requested changes in k-chrispens/sampleworks_figs on
  GitHub origin/main, preserving unrelated work and verifying the exact
  candidate before publication. Invoke only when the user requests landing
  or merging, including an explicit /land invocation, not for review,
  preparation, passing checks, or installing this skill.
disable-model-invocation: true
metadata:
  delta-action: land
---

# Land changes in sampleworks_figs

## Scope and authorization

- Apply only to the repository whose publication remote is
  `github.com:k-chrispens/sampleworks_figs` (SSH or HTTPS).
- An explicit invocation requests completion of the landing workflow.
  Proceed through preparation, verification, publication, and destination
  verification without asking for the same merge permission again.
- The destination is GitHub `origin/main`, not the `local` source-repository
  backlink. Never edit the original checkout or push through `local`.
- Preserve intended existing commits. Commit intended uncommitted changes
  using `<type>(<scope>): <imperative summary>`, at most 72 characters and
  without a trailing period. Do not stage unrelated changes or secrets.
- Resolve conflicts automatically only when the intended result is clear.
  Stop and explain ambiguous conflicts, unresolved scope, failed checks,
  or unsafe changes. Never guess, overwrite unrelated work, force-push,
  rewrite shared history, or bypass protection or review requirements.

## 1. Establish the candidate and current destination requirements

Run from the attached repository's root. Read applicable `AGENT.md`,
`AGENTS.md`, `CLAUDE.md`, contribution policies, and submission templates.
Honor any applicable signing, contributor, human-authorship, review,
documentation, and verification requirements. Do not impose requirements
from the separate Sampleworks repository.

Inspect without discarding changes:

```sh
git --no-optional-locks status --short --branch
git remote -v
git branch -vv
git log -8 --oneline
gh auth status
```

Confirm that `origin` identifies the repository above and that GitHub
authentication permits publication. Do not print credential values.
If intended versus unrelated work is unclear, ask a focused scope question.
Do not silently stash, reset, clean, or stage the whole checkout.

Fetch only after establishing the publication remote:

```sh
git fetch origin
git diff --stat origin/main...HEAD
git diff --name-status origin/main...HEAD
gh api repos/k-chrispens/sampleworks_figs --jq \
  '{default_branch,permissions,allow_merge_commit,allow_squash_merge,allow_rebase_merge,archived}'
gh api repos/k-chrispens/sampleworks_figs/branches/main --jq \
  '{name,protected,commit:.commit.sha}'
gh api repos/k-chrispens/sampleworks_figs/branches/main/protection
gh api repos/k-chrispens/sampleworks_figs/rules/branches/main
gh workflow list --repo k-chrispens/sampleworks_figs --all
```

Interpret errors rather than ignoring them: an authenticated response saying
`Branch not protected`, together with an unprotected branch and no active
rules, is not a verification failure. Authorization, network, or unavailable
rule information is not proof that no requirements exist.

At skill creation this repository used direct pushes to `main`, had no
GitHub workflows, and had no active branch protections/rules. These are
observations, not permanent exemptions. If the default branch or publication
identity changes, stop for clarification. If applicable rules now require
a pull request, use the pull-request path in section 5 rather than a direct
push. Never change rules or use administrative bypass.

## 2. Prepare a clean, isolated landing candidate

Record the starting commit, intended file list, and fetched `origin/main`
commit. Use a uniquely named local branch such as `land/<short-topic>` in
this checkout. Check that the name is unused before creating it.

If there is intended uncommitted work, review its diff, stage only its
explicit paths, check the staged diff, and commit it. Use noninteractive
Git commands; prefix every commit or merge command with `GIT_EDITOR=true`.
The newly installed skill can be included when landing the skill itself
and the originally requested change together.

Before integration, preserve unrelated uncommitted work in place. If it
prevents safe branch creation, integration, or verification, stop and
explain; do not move or discard it without an explicit instruction.

If the fetched destination is not an ancestor of the candidate, merge it
into the isolated landing branch with a normal noninteractive merge:

```sh
GIT_EDITOR=true git merge --no-edit origin/main
```

For a conflict, inspect both sides and their callers/documented intent.
Resolve a clear conflict, stage only the resolution, and complete the merge
noninteractively. If intent is ambiguous, leave the source checkout
untouched, identify the conflicted files, report that nothing has landed,
and wait for the user's decision. Do not drop either side's unrelated work.
Do not use rebase, reset, or force-push to conceal a conflict.

## 3. Verify the exact candidate

Record the candidate commit and check the whole change against the fetched
destination, including any integration resolution:

```sh
git diff --check origin/main...HEAD
git --no-optional-locks status --short
```

Keep verification proportional to the changed files. A documentation-only
change does not require loading large scientific datasets or running
notebooks. Verification results from an earlier candidate do not cover a
subsequent fix, merge resolution, or amended commit: rerun affected checks.

### Documentation, HTML, and skill-only changes

- Read changed prose and ensure links, paths, commands, and examples agree
  with their actual definitions. Do not treat an unimplemented plan's
  example test commands as existing repository tests.
- For changed HTML plans, use a small Python-standard-library
  `html.parser.HTMLParser` check via `uv run --no-project python`, not an
  application dependency installation. Check balanced non-void elements,
  duplicate IDs, internal anchor targets, and required document structure.
  If the document promises self-containment, check that no scripts, styles,
  images, or CSS imports load remote resources. Parsing is not visual
  rendering; state whether the document was opened in a browser.
- For the blueprint HTML, preserve its original-request block and six
  ordered sections unless the requested change explicitly alters them.
  Source: `blueprint/sampleworks-data-contract/plan-sampleworks-data-contract.html`
  defines its header, inline styling, navigation, and section structure.
- For changed skills, parse YAML frontmatter with a real YAML parser when
  available; otherwise stop rather than claiming YAML syntax validation.
  Verify nonempty `name` and `description`, valid referenced paths, and
  for this skill the exact `metadata.delta-action: land` and
  `disable-model-invocation: true` values.
  Runtime source: `uv.lock` includes `pyyaml`; when using that locked
  dependency, use `uv run --locked python` without changing dependency
  declarations or the lockfile.
- Do not run the analysis import smoke check solely for prose/skill changes.

### Analysis-module or dependency changes

Use Python through `uv`; do not invoke a different interpreter or install
with pip/conda. The project requires Python 3.13 or later, with 3.13 selected
by `.python-version`. Use the committed lockfile without updating it:

```sh
uv run --locked python -c "import style, loaders, transforms, plot_scatter, plot_classified, plot_tradeoff, plot_training_bias, plot_budgets"
```

Sources establishing this exact import set and runtime:
`pyproject.toml` defines the Python requirement and analysis dependencies;
`.python-version` selects 3.13; the eight corresponding root `.py` modules
exist and define the imported surfaces; `CLAUDE.md`, Validation, names the
same smoke check. `uv run --help` supports `--locked`.

For changed behavior, also reproduce the affected helper or workflow with
representative inputs and inspect the requested output. If plots change,
consult and update `data_in_each_plot.md` and the notebook figure catalog
as required by `CLAUDE.md`, Adding a New Figure and Reference.

### Notebook changes

Do not execute a marimo notebook as `python <notebook>.py`. The eight
notebooks define `marimo.App` cells, not ordinary script execution.
Sources: their app declarations and `CLAUDE.md`, Commands and Validation.

Verify changed cells through the installed, locked marimo CLI using a
supported export/run command; confirm its exact command and arguments
against that installed version's help before invoking it. Do not copy a
stale documented command without this check. Write any exported executable
script inside the checkout so relative imports resolve; use a unique
temporary filename and remove only that file afterward.

Real-data checks may require ignored `results/` inputs that are absent in
an isolated checkout. Do not bypass their absence by using another thread's
checkout or by modifying the original checkout. If required inputs or the
runtime are unavailable, stop with the missing requirement rather than
claiming the notebook passed. Never run evaluators or patching tools on
the read-only historical source `/mnt/diffuse-shared/raw/sampleworks`.
Sources: `.gitignore` excludes `results/`; `CLAUDE.md` documents the
notebook's data sources and marimo execution model.

No test runner, Makefile, pre-commit configuration, or CI workflow existed
in this repository at skill creation. Reinspect any added test/build
definitions when the change introduces them and run their actual supported
commands; do not invent `pytest`, lint, or build targets.

The interactive-figure refactor adds a standard-library unittest suite and
eight numbered marimo entry points. Run these current checks sequentially:

```sh
uv run --locked python -m unittest discover -s tests -v
uv run --locked marimo check --strict [0-9][0-9]_*.py
```

Sources: `tests/test_*.py` define `unittest.TestCase` tests; `CLAUDE.md`,
Validation, names both invocations; `01_baseline.py` through
`08_sampling_budgets.py` define the current notebook entry points.
Follow `CLAUDE.md`'s remote-runtime requirement: use an isolated pod synced
from the landing checkout, not another thread's checkout or demo workspace.
Confirm source hashes and results-snapshot availability before verification.

## 4. Establish all required remote checks and reviews

Determine requirements from the current branch protection/rules and
applicable repository policies, not solely from which checks happen to
appear. If checks require publishing a candidate branch, publish only a
new, non-destination branch with a normal push, then inspect that exact
candidate's check runs and commit statuses:

```sh
gh api repos/k-chrispens/sampleworks_figs/commits/<candidate-sha>/check-runs
gh api repos/k-chrispens/sampleworks_figs/commits/<candidate-sha>/status
```

Replace `<candidate-sha>` with the actual commit. Account for all pages
when inspecting checks or rules. Wait for every applicable required check
to complete successfully. Pending, failing, cancelled, skipped when
required to execute, missing, or unverifiable required checks block
landing. Verify that any required human review covers the current change;
an automated review is not human approval.

GitHub can return `state: pending` with an empty `statuses` array even when
no status contexts exist. Do not mistake that aggregate for a running
required check. If current rules/policy establish that no remote checks
are required and no CI is configured, report remote verification as
not configured, not passed. Local checks must still pass.

## 5. Publish without bypassing current repository rules

### Normal direct-push path

Immediately before publication, fetch `origin` again. If `origin/main`
advanced, integrate it into the isolated branch as in section 2 and
reverify the resulting candidate and any required checks. Confirm that
the fetched destination is an ancestor of the candidate:

```sh
git merge-base --is-ancestor origin/main HEAD
git push origin HEAD:refs/heads/main
```

A rejected push is not permission to force it. Refetch and safely
integrate a concurrent destination change, or stop if there is a genuine
permission/policy/conflict blocker. Stop if concurrent updates make safe
completion impractical; do not keep retrying indefinitely.

### Pull-request path, only when current requirements demand it

Publish the candidate on a uniquely named branch without force-pushing.
Reuse a matching PR if present. Otherwise supply its branch, title, and
body explicitly so creation cannot open a prompt or editor:

```sh
gh pr create --repo k-chrispens/sampleworks_figs --base main \
  --head "$landing_branch" --title "$pr_title" --body-file "$pr_body_file"
```

Populate those shell variables from the actual prepared branch and the
required submission text; do not send literal placeholders.
Honor submission templates and any human-authored-text requirement:
obtain that text rather than generating it and equating approval with
authorship.

Inspect PR checks, current review requirements, head SHA, and mergeability.
Require every applicable check/review to pass for the actual candidate or
required merge result, including a merge queue when enforced. Use a merge
method permitted by current destination rules/settings, preserving
existing commits where permitted. If requirements force a different
history strategy or human action that cannot be completed safely, explain
the blocker instead of bypassing it.

Perform the actual permitted merge with noninteractive `gh` commands and
wait for the destination update. For a permitted merge-commit strategy,
pin the reviewed head explicitly:

```sh
gh pr merge "$pr_number" --repo k-chrispens/sampleworks_figs \
  --merge --match-head-commit "$candidate_sha"
```

Use a different strategy only when required by the inspected rules and
permitted history policy; never use `--admin`. Merely creating a PR,
enabling a queue, starting checks, or publishing the source branch is
not landing.

## 6. Verify the destination and report accurately

After a successful push or completed PR merge:

```sh
git fetch origin
gh api repos/k-chrispens/sampleworks_figs/branches/main --jq '.commit.sha'
git --no-optional-locks status --short --branch
```

For a direct push, verify the candidate commit is reachable from the
refetched `origin/main` and agrees with GitHub's observed history; a later
concurrent commit is acceptable if it contains the candidate. For a PR
merge, verify the actual merge/squash commit and the requested change at
the destination. Do not require the pre-merge SHA to survive a squash,
and do not claim it did.

Report the target, actual landed commit, checks actually run, remote
check status, and anything not exercised. On a blocker, state that the
changes have not landed and explain the required action. Do not publish
packages, alter branch protections, delete remote branches, or clean up
unrelated files as part of this workflow.
