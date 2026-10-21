---
status: accepted
date: 2026-10-13
decision-makers: Zack (with Claude Code doing the download and storage sessions)
---

# TAHMO gauges kept on one machine; other gauge files on the dev blob with checksums

## Context and Problem Statement

AMMA-CATCH (ADR 0001) only covers three sites in West Africa. To check the products more widely we needed more
gauges that GPCC and CHIRPS don't use. TAHMO's portal offers 26 per-country "corrected and curated 2015-2025
rainfall datasets" (575 stations, mostly East and West Africa) to research users, but its licence forbids sharing.
At the same time, all gauge files had been sitting only in git-ignored `data/` on one laptop, with their provenance
notes, so nothing was reproducible from anywhere else. Two questions: how to take TAHMO, and where gauge data lives.

## Decision Drivers

* Licences: AMMA-CATCH is CC-BY 4.0. TAHMO's research class says *"Research use does not allow for sharing of the
  data to any other person or organization"* and forbids passing data on *"in the form of derived products"* without
  written permission.
* Reproducibility: another person or machine should be able to get the inputs back and confirm they're identical.
* The repo is public: no restricted data, station lists or internal draft correspondence in git.
* Disk on the laptop is limited, but the gauge files are small (AMMA-CATCH 9 MB, TAHMO 640 MB).

## Decision Outcome

1. **TAHMO: all 26 curated country zips, kept exactly as delivered, on one machine only.** They aren't on the team
   blob, Databricks or email. The analysis reduces them locally (`tahmo_daily.py` → `cache/tahmo/`) and publishes
   aggregate scores only (by region and stations-per-cell tier); no station series or coordinates. A checksum list
   is committed (`docs/gauges/tahmo.sha256`), so a re-download from the portal can be verified byte for byte.
   The explainer artifact "Which Rainfall Product Where", which includes TAHMO results, is private for now.
   *Rejected:*
   * Putting TAHMO on the team blob: shared storage counts as sharing under the licence.
   * The portal's "Data export" (per-station exports): the account isn't configured for it, it's capped at about
     30 exports a month, and it gives raw rather than curated data.
   * Asking for the full archive: possible only "after special review", and slower.
   * Downloading only the African files: the non-African ones are tiny and keep the set complete.
2. **Everything else under `data/gauges/` goes on the dev blob** at `projects/ds-precip-intercomparison/raw/gauges/`
   through `scripts/gauges_blob.py`. That's the AMMA-CATCH zips, manifest, download log and notification emails;
   the BWDB station list; and the scouting notes and handover. A sha256 manifest is written last, and `download`
   refuses any file that doesn't match it. The script skips `tahmo/` by name and prints the licence reason, so
   it's never left out silently.
   *Rejected:*
   * Laptop only: a single copy, invisible to anyone else.
   * Committing the files: the repo is public, and data stays out of git.
   * The prod blob: this project writes dev only.
3. **The lasting description lives in git** as `docs/gauges/README.md`: what each source is, licence and citations,
   where the files are, how to get them back, and which scripts use them. It holds no data, no TAHMO station
   details, and none of the internal draft emails; those stay in `SCOUTING.md` on the blob.

### Consequences

* Good: AMMA-CATCH, BWDB and all provenance can be restored with one command, and every file is checked on the way.
* Good: a TAHMO re-download can be proven identical (or shown to have changed) with `shasum -c`.
* Bad: TAHMO exists as one copy on one laptop. Losing it means re-downloading from the portal, which needs an
  active research account (access is time-limited and renewable), and TAHMO may re-issue the files in the meantime.
* Bad: `cache/tahmo/viz.json` holds TAHMO station coordinates for a network map. It must not go into git or onto
  GitHub Pages without TAHMO's permission.
* Open: whether aggregate scores count as "derived products". A question to TAHMO is drafted (in `SCOUTING.md`). If
  they also allow team storage, revisit decision 1 and remove `tahmo` from `LOCAL_ONLY` in `gauges_blob.py`.

### Confirmation

* `python scripts/gauges_blob.py upload`, then `download`, restores 39 files with matching sha256 (run 2026-10-13).
  A second `upload` sends nothing.
* `shasum -a 256 -c docs/gauges/tahmo.sha256` passes on the 26 local TAHMO zips and matches the download log
  written when they were fetched (2026-10-06).

## More Information

* `docs/gauges/README.md`: the map of gauge sources, rules and re-acquisition steps.
* ADR 0001: how AMMA-CATCH was acquired. Its storage decision (git-ignored `data/`) is extended here with the blob copy.
* TAHMO Data Access Policy V1.0: https://usercontent.one/wp/tahmo.org/wp-content/uploads/2024/01/TAHMO-Data-Policy-V1.0.pdf
