# Rain-gauge data used to check the rainfall products

The products in this repo (CHIRPS, IMERG, ERA5, GPCC, …) are compared against each other everywhere, but the
only way to say which one is *right* in a place is to compare them with rain gauges that **none of them used**.
Most national gauges are already inside GPCC and CHIRPS, so we went looking for independent networks. This page
says what we got, where the files are, the rules attached to them, and how to get them back.

The raw files are **not in git** (`data/` is git-ignored, and the repo is public). Decisions and their reasons
are in [ADR 0001](../decisions/0001-amma-catch-gauge-acquisition.md) (AMMA-CATCH) and
[ADR 0002](../decisions/0002-tahmo-gauges-and-gauge-storage.md) (TAHMO, and where gauge data is stored).

## At a glance

| Source | What it is | Size of the network | Period | Licence | Raw files live in | Status |
|---|---|---|---|---|---|---|
| **AMMA-CATCH** | Research rain-gauge networks in Benin, Niger, Mali (plus 3 weather/flux stations incl. Senegal) | 15 datasets; the 10 daily ones hold 355 station records (326 distinct stations) | 1990–2025 (varies by dataset) | CC-BY 4.0 | dev blob `projects/ds-precip-intercomparison/raw/gauges/amma-catch/` | Downloaded 2026-10-05, in use |
| **TAHMO** | Automatic weather stations across Africa, "corrected and curated" rainfall | 575 stations in 26 country files (571 in Africa) | 2015–2025 | TAHMO research licence (no sharing) | **One laptop only**, `data/gauges/tahmo/` | Downloaded 2026-10-06, in use |
| **BWDB** (Bangladesh) | List of the water board's 307 manual rain gauges, with coordinates | 307 stations | 1961–2026 | Public web page | dev blob `…/raw/gauges/bwdb/` | Station list only; rainfall needs a paid request |
| PhilSensors, BWDB data, Ethiopia EMI, Malawi DCCMS, Nepal DHM, FAO SWALIM | National or agency networks | — | — | Request or fee | — | Scouted, not requested. Routes, contacts and draft emails in `SCOUTING.md` on the blob |

## Getting the files onto a machine

```bash
python scripts/gauges_blob.py download   # AMMA-CATCH + BWDB + notes -> data/gauges/, every file sha256-checked
```

`gauges_blob.py upload` sends `data/gauges/` back to the blob (only changed files), with a sha256 manifest written
last. It **skips `tahmo/` on purpose** and prints why.

**TAHMO has to be downloaded again by hand** on any new machine (it can't be shared through the blob, see below):

1. Sign in at https://portal.tahmo.org with a TAHMO research account (Google sign-in). The 2026 download used an
   account that was "not configured", which still sees these files; access under the research class is
   time-limited and renewable.
2. Open **Datasets** ("Corrected and curated 2015-2025 rainfall datasets") and click **Download** on all 26 rows.
   Chrome saves them as `TAHMO_Historical_precipitation_dataset_<Country>.zip`.
3. Put them in `data/gauges/tahmo/` and check they're the same files we used:
   ```bash
   cd data/gauges/tahmo && shasum -a 256 -c ../../../docs/gauges/tahmo.sha256
   ```
   A mismatch means TAHMO has re-issued the files. Results will then differ slightly; note the date.

If the blob copy of AMMA-CATCH were ever lost, ADR 0001 describes the portal steps to extract it again. The portal
keeps updating its data, so a new extraction won't match the old checksums.

## How the analysis uses them

| Script | Reads | Writes |
|---|---|---|
| `scripts/validate_gauges.py` | the 10 daily AMMA-CATCH zips | `site/multivariate/data/gauges.json`, `gauge_viz.json` |
| `scripts/validate_3month.py` | AMMA-CATCH only (via `gauge_viz.json`; borrows helper functions from `validate_tahmo.py`) | `site/multivariate/data/gauges_3month.json` |
| `scripts/tahmo_daily.py` | the 26 TAHMO zips | `cache/tahmo/daily.parquet`, `cache/tahmo/stations.parquet` (local only) |
| `scripts/validate_tahmo.py` | the TAHMO daily cache | `cache/tahmo/scores.json` (aggregate monthly and 3-month scores), `cache/tahmo/viz.json` |

Nothing TAHMO-derived is written under `site/` today; TAHMO results stay in `cache/tahmo/` (git-ignored, local).

`cache/tahmo/viz.json` contains **TAHMO station coordinates** (for a network map). Don't publish it, or anything
built from it, without TAHMO's written permission.

## Rules that come with the data

**AMMA-CATCH**: CC-BY 4.0 (AMMA-CATCH data policy v2.0). Any publication must cite:
1. Galle et al. (2018). AMMA-CATCH a Critical Zone Observatory in West Africa Monitoring a Region in Transition.
   *Vadose Zone Journal*, 17 (1), 18006.
2. AMMA-CATCH (1990). AMMA-CATCH: a hydrological, meteorological and ecological observatory on West Africa. IRD,
   CNRS-INSU, OSUG, OMP, OREME. doi:10.17178/AMMA-CATCH.all

Dataset DOIs are optional extras (listed in the AMMA-CATCH `MANIFEST.md` on the blob). Offer co-authorship to, or
at least acknowledge, the data providers.

**TAHMO**: research class of the TAHMO Data Access Policy V1.0
(https://usercontent.one/wp/tahmo.org/wp-content/uploads/2024/01/TAHMO-Data-Policy-V1.0.pdf). In its own words:
- *"Research use does not allow for sharing of the data to any other person or organization."*
- *"Data accessed from TAHMO cannot be sold or passed on to other users without the explicit permission from TAHMO
  in writing whether in their original or modified formats or in the form of derived products."*
- Cite: Van de Giesen, N., Hut, R. and Selker, J. (2014). The trans-African hydro-meteorological observatory
  (TAHMO). *WIREs Water* 1(4): 341–348. doi:10.1002/wat2.1034.

In practice:
- Keep the files on one machine, and off the team blob, Databricks and email.
- Publish only aggregate scores (by region or station tier), never station series or coordinates.
- Whether even aggregate scores count as a "derived product" is still to be confirmed with TAHMO (see open items).

## TAHMO in more detail

Each country zip holds `metadata.csv` (station ID, name, country, latitude, longitude, elevation, gauge height,
installation date) and one CSV per station with `time, station_id, precip_mm, precip_sensor_id,
final_quality_flag`. Points to know:
- **Time step:** 5-minute rainfall (10 min in Ethiopia and Madagascar, 15 in Niger).
- **Time zone:** not stated; `tahmo_daily.py` treats it as UTC.
- **Missing values:** about 1% of rows are empty or `nan`.
- **Quality flags:** every row has flag 1.0, so TAHMO appears to have removed flagged data already, and a gap means
  missing data. The flag scheme isn't documented in the files.
- **Overlap with the analysis period:** most stations were installed in 2015–2019, so overlap with 2001–2020 is
  only about 2015–2020; the TAHMO scores use 2015–2025.

| File suffix | Region | Stations | Records span | Step (min) |
|---|---|---|---|---|
| Kenya | East | 116 | 2015-03 → 2025-09 | 5 |
| Uganda | East | 48 | 2015-07 → 2025-09 | 5 |
| Tanzania | East | 24 | 2017-03 → 2025-09 | 5 |
| Rwanda | East | 14 | 2015-11 → 2025-09 | 5 |
| Ethiopia | East | 6 | 2022-06 → 2025-09 | 10 |
| Malawi | Southern | 23 | 2017-12 → 2025-09 | 5 |
| Zambia | Southern | 12 | 2018-07 → 2025-09 | 5 |
| Mozambique | Southern | 11 | 2018-07 → 2025-09 | 5 |
| SouthAfrica | Southern | 7 | 2017-07 → 2025-09 | 5 |
| Madagascar | Southern | 1 | 2023-07 → 2025-09 | 10 |
| Zimbabwe | Southern | 1 | 2018-09 → 2019-04 | 5 |
| Ghana | West | 90 | 2015-04 → 2025-09 | 5 |
| Nigeria | West | 41 | 2016-06 → 2025-09 | 5 |
| Mali | West | 30 | 2016-04 → 2025-09 | 5 |
| Togo | West | 18 | 2017-09 → 2025-09 | 5 |
| Benin | West | 10 | 2017-03 → 2025-09 | 5 |
| BurkinaFaso | West | 9 | 2017-04 → 2025-09 | 5 |
| Niger | West | 9 | 2019-06 → 2025-09 | 15 |
| Cote_D_Ivoire | West | 4 | 2019-11 → 2025-09 | 5 |
| Senegal | West | 4 | 2022-09 → 2025-09 | 5 |
| CDR | Central | 8 | 2018-01 → 2025-09 | 5 |
| Cameroon | Central | 7 | 2018-09 → 2025-09 | 5 |
| Other | Mixed (21 countries) | 78 | 2015-03 → 2025-09 | 5 |
| USA | Outside Africa | 2 | 2015-09 → 2025-09 | 5 |
| CostaRica | Outside Africa | 1 | 2019-11 → 2024-06 | 5 |
| Germany | Outside Africa | 1 | 2019-12 → 2025-09 | 5 |

**Checks done on download (2026-10-06):**
- Every zip passed its integrity test, and its station list matches its station files exactly.
- Country codes match the file names, except `Other`, which is mixed by design.
- No station appears in two files, and every station file has the expected columns.
- The sha256 checksums are in [`tahmo.sha256`](tahmo.sha256).

Twenty files downloaded normally. Six (Other, Togo, Uganda, USA, Zambia, Zimbabwe) had to be re-requested. Those six
were recovered from Chrome's temporary files, identified by the country codes inside them, and copied byte for
byte.

## AMMA-CATCH in more detail

The full record is in `MANIFEST.md`, `download_log.tsv` and `emails/` inside the blob folder, and the reasoning is
in ADR 0001. What matters most when using it:

- **Daily totals cover the 24 hours ending at 06:00 UTC**, so a value stamped day D is mostly rain from day D−1.
- **Negative values mean missing.** A negative value other than −9999.9 is the total for a whole gap, estimated
  from a collection bucket, and is repeated on each day of the gap.
- **The two Niamey datasets overlap.** `CE.Rain_Nc` and `CE.Rain_Nc_9004` share stations, and each zip contains a
  few 2004/2005 files labelled with the other dataset's code.
- **A few stations are missing or ambiguous.** Two `CL.Rain_GT` stations share the name GOURMA_RHAROUS. `AL.Met_G`
  has no rainfall for KOBOU, and `AE.H2OFlux_Ncw` has rainfall for only 1 of its 4 stations.
- **A few coordinates change between years.** Cotonou's 2016–2018 files place the station about 175 km offshore.
- **Every download notified the dataset's PI.** The portal asks users to propose collaboration.

## Open items

- **Write to TAHMO** (info@tahmo.org) to confirm two things: that publishing aggregate validation scores is allowed,
  and whether the files may sit on the team's private storage. A draft is in `SCOUTING.md` on the blob. Until they
  reply, TAHMO stays on one machine and nothing station-level is published.
- **TAHMO access lapses** unless it's renewed; renewal needs evidence of research use.
- **Request-only sources** (PhilSensors is the strongest: about 1,400 rain stations across the Philippines) are
  waiting on requests. Routes, costs and draft emails are in `SCOUTING.md` on the blob.
