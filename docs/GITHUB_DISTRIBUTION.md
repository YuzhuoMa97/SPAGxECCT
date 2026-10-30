# Source and GitHub Releases

| User | Download | Installation |
|---|---|---|
| Historical R workflow | R source with compatible legacy dependencies | Existing R interface |
| C++ source user | Clone, Code ZIP or tagged source archive | Read root README_CPP.md, build inside cpp/ |
| Windows binary user | Complete Windows ZIP from Releases | Extract, run cpp/SPAGxE.exe or root run_examples.cmd |

R installation does not install the independent EXE. R build exclusions omit
cpp/docs from the R tarball, not the GitHub repository/source ZIP.

A commit records source changes. A tag names a commit. A Release supplies a
version page, notes and assets for a tag. GitHub automatically offers tagged
source ZIP/tar; maintainers separately attach compiled software. A source push
is enough for source availability, but not for the precompiled download route.
A Release is not a second algorithm mode.

## Validate develop before merging main

1. Keep existing history/R files and add complete cpp/, docs/, notices and CI.
   Preserve the clone's own .git metadata. Save existing personal changes and
   fetch/pull the intended base before applying a prepared source update.
2. In GitHub Desktop, create or select develop from an up-to-date main. Check
   Current Branch before copying files. Review and commit the complete addition,
   then Publish branch or Push origin. Include vendored dependencies and
   examples, not only cpp/src/.
3. On GitHub, select develop and record its exact commit SHA. Check the Actions
   results for that commit and each intended platform. The workflow accepts
   pushes on any branch, including develop, and pull requests; relevant source,
   documentation and workflow paths trigger the checks.
4. Select develop before Code -> Download ZIP, or download the recorded commit
   archive for a fixed validation target. Extract into a fresh directory,
   build inside cpp/ using the installation instructions, and run the version
   check and all six tutorials without reusing previous build objects.
5. If an Actions candidate bundle is used, verify its commit record and archive
   checksum, extract the complete package, and run its root example launcher.
   Record the downloaded source commit, platform, commands and outcomes.
6. Open a pull request with base=main and compare=develop. Review the diff and
   required checks, resolve conflicts, then merge. Check the resulting main
   commit and its Actions results. If merging changes source, dependencies,
   build settings or examples, rebuild and revalidate the corresponding bundle.

In a public repository, develop is also public. It isolates development changes
from main; it is not storage for confidential code or data. Private data, caches,
internal audits and unpublished-method implementation materials stay outside
the public addition.

## Candidate artifacts and platform evidence

The workflow builds, runs numerical/reader regressions, runs the six tutorials,
creates a complete bundle with make dist, extracts that archive, and reruns the
bundle tutorials. It checks all 45 non-log tutorial outputs against the source
run. Windows bundle examples are invoked through native cmd with a system-only
PATH. Each successful job uploads candidate-SPAGxE-* with the commit identifier,
the complete tar.gz bundle and its SHA256 record. The Actions artifact download
is a wrapper: extract it first, verify the enclosed tar.gz, then extract that
complete bundle and run run_examples.cmd or run_examples.sh for its platform.

These are commit-specific test candidates. They are not automatically published
Releases and do not replace the maintainer's verified Windows Release ZIP.
The separate validation-* artifacts contain job evidence. A workflow file is
not proof that a platform run passed; inspect the actual job for the target
commit before making platform claims.

The initial develop push can trigger CI directly. Manual Run workflow becomes
available only after the workflow exists on the default branch; when available,
select the intended branch explicitly. Documentation path filtering also matters:
an absent or skipped run is not a successful validation result.

## Create the matching Release

1. Open Releases -> Draft a new release after the final main commit is accepted.
2. Create cpp-v1.2.0 targeting the exact commit containing the accepted source,
   and describe the methods, platform requirements and verification limits.
   Retain the historical R tag. If this C++ tag was already published, use a
   consistently updated software version and tag instead of replacing its code
   or assets with a changed algorithm.
3. Attach SPAGxE-1.2.0-windows-x86_64.zip and SHA256SUMS.txt listing that exact
   asset. The complete bundle includes the executable, examples, attribution
   and licenses. Confirm it corresponds to the accepted source and build policy.
4. Use the tag's automatic Source code ZIP/tar as the corresponding source;
   do not substitute an unrelated older snapshot. A source archive contains
   build inputs, while the prepared Windows ZIP contains the runnable program.
5. Save a draft and review the target, tag, English notes and both attachments.
   Initially mark pre-release while broader platform evidence remains pending.
   If immutable releases are enabled, attach all assets before publication.
6. Publish only when the claims match the evidence. Download the actual assets,
   verify SHA256, extract the complete folder, and run the version/tutorial
   checks. Also confirm the tagged source contains complete cpp/ and builds as
   documented. Record the final commit and asset checksums with the release.

CI has read-only repository permissions and does not create tags, publish
Releases, or upload the maintainer's final Release assets automatically.

Official instructions:
[Desktop commit/push](https://docs.github.com/en/desktop/making-changes-in-a-branch/committing-and-reviewing-changes-to-your-project-in-github-desktop),
[Desktop branches](https://docs.github.com/en/desktop/making-changes-in-a-branch/managing-branches-in-github-desktop),
[Desktop pull requests](https://docs.github.com/en/desktop/working-with-your-remote-repository-on-github-or-github-enterprise/creating-an-issue-or-pull-request-from-github-desktop),
[Source archives](https://docs.github.com/en/repositories/working-with-files/using-files/downloading-source-code-archives),
[Workflow events](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows),
[Release creation/assets](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository).

See [installation](INSTALLATION.md), [release scope](RELEASE_STATUS.md) and
[ownership/attribution](../OWNERSHIP.md).
