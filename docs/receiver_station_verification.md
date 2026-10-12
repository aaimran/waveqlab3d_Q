# Receiver station implementation verification

Verified October 11, 2026, using the GNU Fortran Debug build with bounds checks
and floating-point traps enabled. The implementation supports the separate
`&output_receiver_stations` namelist, explicit `append_block`, and per-station
`interface_stations` selection. See [configuration guide](receiver_station_output.md).

## Checks completed

- Debug build: passed.
- `receiver_station_output`: all 54 cases passed (213.09 seconds), including named stations.
- Existing `cQ_cgt_traditional6_dynamic`, `cQ_cgt_traditional6_histories`,
  `q8_dynamic_decomposition`, and `q4_dynamic_decomposition`: passed.
- Legacy station files compared with the unchanged executable from repository
  HEAD: all four files matched byte for byte in both one- and two-rank runs.
- Full-fault outputs compared with that unchanged executable: all eight files
  matched byte for byte for both Cartesian and curved faults.
- Source diff whitespace check: passed.

## Receiver coverage

The suite covers new and legacy parsing; both interface sides and individual
side selection; off-interface stations in both blocks; nearest-sample writer
selection across 1–4 ranks; numbered, coordinate and grid-index names; separate
block folders; collisions; external station lists; empty and disabled output;
headers and metadata; column reordering; case normalization; curved interfaces;
malformed rows and sections; conflicting controls; missing external lists;
unmapped stations; unavailable block 2; duplicate IDs; directory failures; and
paths containing spaces. A slipping case proves that side selection returns
different velocity histories. A locked CG-T case checks numerical agreement
across MPI layouts. Named and numbered histories match at the same MPI layout. Named tests
also verify the four actual TPV5 filenames, external lists, metadata, case-insensitive
duplicate detection, invalid and incomplete rows, suffixes and exact-file
collisions. New slipping histories match legacy histories at each
layout, rather than requiring the existing slipping solver to be invariant
across layouts.

## Existing limitations

The older Cartesian and curved full-fault truth fixtures did not match the
current executable. Their default station paths also refer to missing nested
folders. Disabling station output let their full-fault comparisons run, but
the slip differences against old truth files were 1.9621861 m and 5.847215 m,
respectively. The unchanged executable produces byte-identical full-fault
files, establishing that these differences predate this implementation.

`output_exact_moment` still creates the selected files and headers, but the
solver's analytic history calls are commented out. Its filename/header
behavior was verified; analytic samples are not implemented by this change.
Nearest-node receiver mapping and TPV5 displacement/upload conversion retain
their existing behavior.

## Ready-to-run build

The verified executable is `bin/waveqlab3d-receivers-debug`. Generated changes
in the tracked `build-cgt-debug` directory are excluded from the source patch;
reconfigure and rebuild that directory to regenerate its new test registration.
The executable SHA-256 is `32e6f57215289db0f6a334418421060c9c09dd8f787bda813d9cba8a91a131cc`.
