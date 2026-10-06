---
status: accepted
date: 2026-10-05
decision-makers: Zack (with Claude Code doing the download session)
---

# Acquire AMMA-CATCH rain gauges as daily, verified, renamed zips kept git-ignored

## Context and Problem Statement

Validating the gridded products needs rain gauges that are dense and, ideally, not already in GPCC/CHIRPS.
AMMA-CATCH (Benin, Niger, Mali, Senegal) fits, but its portal (http://bd.amma-catch.org) has no API: you
log in, build one extraction per dataset (stations × variable × period × options), accept the data policy,
and the archive arrives later as an emailed link that only works with a logged-in session, expires after a
week, and is served over plain http. Several choices had to be made: what to fetch, at which time step, how
to get the files onto disk, and how to name, verify and store them so another session could use them
unattended. The result is `data/gauges/amma-catch/` (see its `MANIFEST.md`), read by
`scripts/validate_gauges.py`.

## Decision Drivers

* Complete, documented coverage: a validation set that is silently missing stations or years is worse than useless.
* Fail loudly: every file is checked against what was requested; discrepancies are recorded, not smoothed over.
* Licence CC-BY 4.0 with mandatory citations, and PI notification on every download: keep provenance with the data.
* Data stays out of git (`data/` is ignored); the consumer is another session that finds files by name.
* No credentials handled by the agent; licence acceptance is the user's act (or explicitly delegated).

## Decision Outcome

Eleven linked choices. For each: what was chosen, then the rejected alternatives and why.

1. **Scope: all 15 datasets in the portal's *Precipitation* category.** That is 12 rain-gauge datasets (core)
   plus 3 flux/meteo stations (`AE.H2OFlux_Ncw`, `AL.Met_G`, `PA.Met_SNNs`) taken for their precipitation
   variable only, marked secondary.
   *Rejected:* only `CL.Rain_*` codes, which would miss the `CE.*` Niamey sets (`CE.Rain_Nc_9004` alone has
   166 stations) and the daily-only `*.RainD_*` sets. Also rejected: all variables of the meteo datasets,
   which are irrelevant here.

2. **Time step: the provider's daily variable** ("Precipitation Amount (previous 24 hours)", 24 h ending
   06:00 UTC) wherever it is offered. The daily-only datasets are taken as they are, and the station
   datasets at their native step (15 min, 30 min, or an mm/h rate).
   *Rejected:* the 5-minute variable, which is much larger and would have to be aggregated here,
   re-implementing the provider's gap and bucket handling. Also rejected: hourly, which brings no benefit
   for monthly validation.

3. **Stations and period: every station explicitly ticked, full period.** The portal's default selection is
   partial for several datasets (only 1–2 of Donga's 50 stations and of `CL.Rain_N`'s 31 were ticked). Its end-date validator
   caps the period at the last data date.
   *Rejected:* accepting the defaults, which would have dropped most stations without any error.

4. **Extraction options: all four ticked** (quality flag, computation flag, sensor info, measurement height).
   The quality flag is what separates `checked data` from `missing value` and from
   `Cumulative rainfall estimated from control bucket` (negative gap totals).
   *Rejected:* the defaults (no flags), which leave estimates and real values indistinguishable.

5. **Format: "CSV: 1 file by station and by year".**
   *Rejected:* one CSV per variable-year with all stations (wide, and harder to check per station), NetCDF
   (CSV/text was asked for, and the per-file headers carry station metadata and licence), and ORE (annual
   totals only).

6. **Delivery: through the user's logged-in Chrome session.** Chrome's download folder was pointed at
   `data/`. Each zip was copied from Chrome's own file once its zip CRC passed. Often that was still the
   held temp file (`Unconfirmed *.crdownload`, `.com.google.Chrome.*`) behind Chrome's insecure-download or
   Save-As prompt, because the bytes are complete before the prompt is answered. The user allowed
   "Automatic downloads" and "Insecure content" for the site, without which Chrome blocked every download
   after the first.
   *Rejected:*
   * `curl` with the session cookie: the emailed links need the authenticated session, and pulling the
     session token out of the browser was not done.
   * A localhost receiver fed from the page: expected to be blocked by Chrome's local-network-access rules for public http pages (not tried).
   * Dragging files from `~/Downloads`: macOS privacy blocks the terminal from reading it, and it costs the
     user more manual steps.

7. **Licence acceptance: done by the user, then explicitly delegated.** The user ticked the first boxes and
   then told the agent to tick the rest. The policy (v2.0, June 2025, CC-BY 4.0) and its mandatory citations
   are reproduced in `MANIFEST.md`.

8. **File names: renamed to `AMMA-CATCH_<code>_<tag>.zip`.** The tag is `daily` for the 24 h variable;
   otherwise it is the variable name, e.g. `_Precipitation_Amount`. The portal's own name is
   `<code>_BD_AMMA-CATCH_<YYYY_MM_DD>.zip`.
   *Why:* the portal name carries the extraction date but not the variable, so a re-download gives a
   different name and a second extraction of the same dataset (e.g. 5-minute) would be indistinguishable.
   Most files were also captured from nameless Chrome temp files.
   *Rejected:* keeping portal names (date-based, no variable, and unknown for temp captures) and unzipping
   into a CSV tree (the user asked to keep the zips; about 3,900 small files).
   *Cost:* see Consequences, since a code is a prefix of another code.

9. **Verification and logging: each zip checked against its notification email.** The checks are dataset
   code, variable, station list, period, CRC, and a sha256 recorded in `download_log.tsv`. These raise an
   error: a corrupt zip, a member file with an unexpected dataset code, or a variable mismatch. Merely
   recorded (column `issues`): stations missing or extra, and years outside the requested period. Members
   labelled with another dataset's code are accepted only when that code is named explicitly, and are then
   listed by name. Emails are saved in `emails/`, with the Outlook safelinks wrapper replaced by the portal
   URL. The `CE.Rain_Nc_9004` email was pasted truncated, so its expected values came from the portal panel
   the request was built from (`emails/CE.Rain_Nc_9004.expected.txt`).
   *Rejected:* trusting file names or the portal's success message. The first batch of launches looked fine
   but silently did nothing, most likely because of stale per-tab JSF view state (reloading fixed it), and only the per-tab confirmation
   message caught it. Also rejected: dropping cross-labelled members quietly.

10. **Storage and hand-off: `data/gauges/amma-catch/`, git-ignored.** It holds the zips (not unzipped), the
    `MANIFEST.md` (what, where, licence, citations, gaps), `download_log.tsv` and `emails/`, plus an empty
    `DONE` file created last as the signal to the analysing session.
    *Rejected:* committing data (licence attribution and repo size) and leaving files in `data/` under
    Chrome's names.

11. **Cleanup: Chrome's own saved copies in `data/` were deleted only after their sha256 matched the logged
    copy.** Held temp files were left for the user to discard in Chrome, since deleting them under an open
    prompt was not worth the risk.

### Consequences

* Good: all 15 extractions arrived and passed verification (about 9 MB in total), each with a provenance
  trail: email, sha256, and expected-versus-found counts.
* Good: the daily convention is explicit. A value stamped D 06:00 UTC is the 24 h ending then, and
  `validate_gauges.py` assigns it to day D−1.
* Bad: **name collision in a downstream glob.** `validate_gauges.py` finds zips with
  `glob(f"AMMA-CATCH_{code}_*.zip")` and takes `next()`. For `CE.Rain_Nc` that also matches
  `AMMA-CATCH_CE.Rain_Nc_9004_daily.zip`. It picks the right file today only because of directory listing
  order. If the order flips, the 2004–2025 `CE.Rain_Nc` series silently disappears. **Fixed 2026-10-05:**
  `load_stations()` now requires an exact name and raises unless exactly one zip matches (results unchanged).
* Bad: portal metadata quirks are carried through as delivered (details in `MANIFEST.md`):
  * `CE.Rain_Nc` and `CE.Rain_Nc_9004` each contain a few 2004/2005 files labelled with the other code.
    Same-named files differ between the zips (KALASSI-2004 by 0.01 mm on one day, and the 2005 files in
    `_9004` are single −9999.9 rows).
  * Two distinct `CL.Rain_GT` stations share the name `GOURMA_RHAROUS`, so one may be hidden in years
    where both reported.
  * `AL.Met_G` KOBOU and three of the four `AE.H2OFlux_Ncw` stations have no precipitation files.
* Bad: station coordinates change between yearly files, and consumers that keep the first coordinates they
  read depend on member order:
  * `CL.Rain_Cotonou` (ASECNA) puts the station 175 km offshore (5.412N 1.116E) in 2016–2018, against
    6.357N 2.388E in 2008–2015.
  * Eleven `CL.Rain_N` stations shift 0.2–4.3 km from 2001 (GUILAHEL from 1995). None changes 0.5° cell,
    but KOURE_SUD (4.3 km), BERKIAWEL (3.7), DAREY (3.6) and GUILAHEL (2.5) move further than the 2 km CHIRPS station-matching threshold.
* Bad: re-acquiring means repeating the portal workflow (the emailed links expire after a week). Every
  download emails the dataset PI, and the portal reminds users to propose a collaboration.
* Neutral: the station datasets use other quality-flag wording (`validated data`, `unchecked data`,
  `Correct value (all checks passed)`). `validate_gauges.py` excludes them now and would raise on the
  unknown flags if they were added.

### Confirmation

* `download_log.tsv` has exactly 15 rows, one per dataset code. A re-run of the sha256 check against the
  zips on disk passed after the last file change.
* The check script and its runs lived in the session's scratchpad, not the repo. The log, emails and
  manifest are the durable record.

## More Information

* `data/gauges/amma-catch/MANIFEST.md`: per-dataset stations, period, time step, DOI, citation, PI, notes and gaps.
* AMMA-CATCH data policy: http://bd.amma-catch.org/resources/data_policy/data_policy_BD_AMMA-CATCH_en.pdf
* Revisit if the 5-minute data are needed (sub-daily validation): extract it as a separate variable. The
  naming in decision 8 keeps it distinct (`AMMA-CATCH_<code>_Precipitation_Amount_previous_5_minutes.zip`).
