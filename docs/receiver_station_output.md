# Receiver station output

Use `&output_list` for exact-moment and full-field controls, and
`&output_receiver_stations` for receiver configuration. Existing input files
with receiver settings in `&output_list` remain supported. Do not supply
receiver controls in both sections: preflight rejects that ambiguity, even if
both values agree. Exactly one `&output_list` is required and the receiver
section is optional. Repeated sections and malformed sections are errors.

## TPV5 off-fault example

```fortran
&output_list
 output_exact_moment = F,
 output_fault_topo = F,
 output_fields_block1 = F,
 output_fields_block2 = F,
 stride_fields = 1
/

&output_receiver_stations
 output_seismograms = T,
 output_station_info = F,
 output_station_mapping = F,
 station_xyz_index = F,
 station_list = 'infile',
 station_list_file = '',
 station_file_directory = 'output/tpv5_off_fault',
 station_output_order = 't vx vy vz',
 station_number_in_list = T,
 station_number_in_filename = T,
 station_use_block_subdirectories = F,
 append_block = F,
 interface_stations = 'block_2',
 station_add_header = T,
 station_add_metadata = T
/

!---begin:station_list---
1  3.0  0.0  -12.0
2  3.0  0.0   12.0
3  3.0  7.5  -12.0
4  3.0  7.5   12.0
!---end:station_list---
```

For the repository's TPV5 geometry, coordinates are `(fault-normal, depth,
strike)` in km. With a problem name of `tpv5`, the numbered filenames are
`output/tpv5_off_fault/tpv5_station-1.dat` through `tpv5_station-4.dat`.
The four numerical columns are time and Cartesian velocities. This does not
add displacement integration or SCEC upload conversion.

## Named station rows

Set `station_name_in_list=T` in `&output_receiver_stations` to use
`name x y z` rows. Set `station_number_in_list=F`,
`station_number_in_filename=F`, and `station_xyz_index=F` in this mode.
Named filenames use the supplied name without the problem-name prefix.
For TPV5:

```fortran
 station_name_in_list = T,
 station_number_in_list = F,
 station_number_in_filename = F,
 station_xyz_index = F,
 append_block = F,
 station_use_block_subdirectories = F,
```

```text
!---begin:station_list---
body030st-120dp000  3.0  0.0  -12.0
body030st120dp000   3.0  0.0   12.0
body030st-120dp075  3.0  7.5  -12.0
body030st120dp075   3.0  7.5   12.0
!---end:station_list---
```

This writes `body030st-120dp000.dat`, `body030st120dp000.dat`,
`body030st-120dp075.dat`, and `body030st120dp075.dat` in the configured
station directory. Names must have 1–64 ASCII letters, digits, hyphens or
underscores; duplicate names are rejected without regard to letter case.
Names with directory separators or dots are rejected. The default is F,
so existing numbered and coordinate-only lists keep their behavior.
External station lists support the same named-row format. Metadata includes
`# station_name: <name>` when enabled.

`append_block=T` adds `_block1` or `_block2` before `.dat`, and block
subdirectories remain independent. Exact-moment filenames are
`<station_name>_exact[_blockN].dat`; collisions with any other selected
station filename are rejected before files open.

## Interface selection and filenames

`interface_stations` accepts `block_1`, `block_2`, or `both` (default), ignoring
case and surrounding whitespace. It applies to stations whose nearest sampled
node is on the shared q face: q=nq in block 1 and q=1 in block 2. This uses mesh
indices, including for curved interfaces. Mapping remains nearest-node sampling
with the existing search radius; it does not add interpolation or exact
geometric projection. Check the mapped coordinate and distance in metadata.
Interior stations remain eligible in their own block, regardless of this setting.
`block_2` requires a two-block model.

`append_block` defaults to T for the new format and appends `_block1` or
`_block2` immediately before `.dat`, in all three naming modes:

- Number: `<name>_station-<number>[_blockN].dat`
- Coordinate: `<name>_<x>_<y>_<z>[_blockN].dat` (three decimal places)
- Grid index: `<name>_<i>_<j>_<k>[_blockN].dat`

`station_use_block_subdirectories` independently controls `block1/` and
`block2/` directories. Both interface sides can be written with
`append_block=F` if these subdirectories keep paths distinct. Duplicate paths
across stations, blocks, or MPI ranks are rejected before station files open.
New-format stations with no selected mapped node are also rejected. Each selected
block uses one MPI writer per station, choosing the nearest eligible sample and
then the lowest Cartesian communicator rank for an exact tie.

The remaining defaults match the legacy receiver controls: seismograms disabled,
station information and mapping messages enabled, internal station list, output
folder `seismogram`, column order `t vx vy vz`, grid-index naming, no station
numbers in rows, block subdirectories enabled, headers and metadata disabled.
`station_xyz_index` chooses coordinate filenames, not grid-index input rows.

Legacy `common_stations_blocks='block1'|'block2'|'both'` remains supported only
in the old format and retains its block-wide selection and automatic suffix
behavior. It is not an alias for `interface_stations`. Exact-moment receiver
files share the station configuration and filename rules and retain their
`_exact` filename marker. The existing analytic exact-moment time-history calls
are commented out in the solver; enabling this option currently creates files
and headers, without analytic samples. This change preserves that behavior.
Full-field output and fault-receiver output keep their
existing controls.

## Verification

`tests/receivers/check_outputs.py` runs small slipping- and locked-interface
cases. It checks side selection, actual sample values against legacy output,
MPI layouts, filename modes, external lists, unmapped stations, and rejection
of collisions and malformed configuration. Run it through the
`receiver_station_output` CTest or directly with `--exe`, `--mpiexec`, and
`--work`. Use a fresh work directory for each run.

The slipping regression compares new and legacy histories at each MPI layout:
the existing slipping solver has serial/MPI numerical differences. The locked
CG-T case independently checks agreement across 1, 2, 3 and 4 ranks. Receiver
configuration changes do not address the slipping solver's numerical differences.
